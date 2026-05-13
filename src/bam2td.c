/*
 *   bam2td: HOMER tag directory creator from BAM/SAM/CRAM via htslib.
 *   Copyright (C) 2026  Benjamin Jean-Marie Tremblay
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 */

// bam2td reads a coordinate-sorted BAM/SAM/CRAM via htslib and writes a HOMER
// tag directory consumable by downstream HOMER tools (findPeaks, annotatePeaks,
// analyzeRepeats, ...). The output is:
//
//   <out>/<chr>.tags.tsv    5' tag records, one line per (pos, strand) bucket
//   <out>/petag.tsv         fragment records, paired-end only
//   <out>/tagInfo.txt       summary header consumed by HOMER
//
// The implementation is single-pass and never holds reads in RAM. Coordinate-
// sorted input is required so adjacent identical (pos, strand) entries can be
// collapsed on the fly.

#include <errno.h>
#include <inttypes.h>
#include <libgen.h>
#include <locale.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "version.h"

#define DEFAULT_MAPQ          10
#define DEFAULT_MIN_LEN        0
#define DEFAULT_MAX_LEN        0   // 0 disables the cap
#define DEFAULT_TBP            0   // 0 disables the cap
#define DEFAULT_PRECISION      1
#define DEFAULT_MIS           -1   // -1 disables NM filtering
#define DEFAULT_THREADS        1
#define DEFAULT_FRAGLEN      150

// Diagnostics ---------------------------------------------------------------

#define quit(...) do {                                       \
    fprintf(stderr, "[E::%s] ", __func__);                   \
    fprintf(stderr, __VA_ARGS__);                            \
    fputc('\n', stderr);                                     \
    exit(EXIT_FAILURE);                                      \
  } while (0)

#define warn(...) do {                                       \
    fprintf(stderr, "[W::%s] ", __func__);                   \
    fprintf(stderr, __VA_ARGS__);                            \
    fputc('\n', stderr);                                     \
  } while (0)

#define msg(...) do {                                        \
    if (opts->v) fprintf(stderr, __VA_ARGS__);               \
  } while (0)

static void *alloc(size_t size) {
  void *p = calloc(1, size);
  if (p == NULL) quit("Out of memory (requested %zu B).", size);
  return p;
}

static long get_mem(void) {
  struct rusage r;
  getrusage(RUSAGE_SELF, &r);
  long bytes = r.ru_maxrss;
#ifdef __linux__
  bytes *= 1024;                  // ru_maxrss is KB on linux, bytes on macOS
#endif
  return bytes;
}

static void print_mem(FILE *out) {
  fprintf(out, "Peak memory usage: %'.2f MB\n",
          ((double) get_mem() / 1024.0) / 1024.0);
}

static void print_time(FILE *out, time_t s) {
  if (s > 7200) {
    fprintf(out, "%'.2f hours.\n", ((double) s / 60.0) / 60.0);
  } else if (s > 120) {
    fprintf(out, "%'.2f minutes.\n", (double) s / 60.0);
  } else {
    fprintf(out, "%'ld second%s.\n", s, s == 1 ? "" : "s");
  }
}

// Options -------------------------------------------------------------------

typedef struct opts_t {
  const char *input;
  const char *output_dir;
  const char *genome;
  int  mapq;
  int  min_len;
  int  max_len;
  int  tbp;
  int  precision;
  int  mis;
  int  threads;
  int  default_fraglen;
  bool flip;
  bool single;
  bool sspe;
  bool rmsoft;
  bool omitSN;
  bool keep_all;
  bool keep_one;
  bool unique;
  bool read1;
  bool read2;
  bool keep_dup;
  bool keep_qcfail;
  bool v;
} opts_t;

static opts_t default_opts(void) {
  opts_t o = {0};
  o.mapq            = DEFAULT_MAPQ;
  o.min_len         = DEFAULT_MIN_LEN;
  o.max_len         = DEFAULT_MAX_LEN;
  o.tbp             = DEFAULT_TBP;
  o.precision       = DEFAULT_PRECISION;
  o.mis             = DEFAULT_MIS;
  o.threads         = DEFAULT_THREADS;
  o.default_fraglen = DEFAULT_FRAGLEN;
  return o;
}

static void usage(FILE *out) {
  fprintf(out,
"bam2td " BAM2TD_VERSION " - HOMER tag directory from BAM/SAM/CRAM\n"
"\n"
"Usage:\n"
"  bam2td <output-dir> <input.bam|sam|cram> [options]\n"
"\n"
"Options:\n"
"  -genome <name>       genome name written into tagInfo.txt (default unknown)\n"
"  -mapq <int>          minimum MAPQ (default %d)\n"
"  -minlen <int>        minimum tag/fragment length (default %d)\n"
"  -maxlen <int>        maximum tag/fragment length (0 = no cap)\n"
"  -tbp <int>           cap reads per (pos,strand) (0 = no cap)\n"
"  -precision <1|2|3>   decimal places when printing counts (default %d)\n"
"  -mis <int>           maximum NM:i (mismatches+indels); -1 disables\n"
"  -flip                flip strand of all reads\n"
"  -sspe                strand-specific PE: flip strand of R2\n"
"  -single              treat input as single-end (skip petag.tsv)\n"
"  -rmsoft              report query length excluding soft-clipped bases\n"
"  -omitSN              omit unused @SQ entries from tagInfo.txt\n"
"  -keepAll             keep secondary and supplementary alignments\n"
"  -keepOne             keep primary alignments only (the default)\n"
"  -unique              require NH:i==1 where present\n"
"  -read1 / -read2      in PE data, keep only R1 or R2\n"
"  -keepDup             do not skip BAM_FDUP reads (default: skip)\n"
"  -keepQCfail          do not skip BAM_FQCFAIL reads (default: skip)\n"
"  -p <int>             htslib decoder threads (default %d)\n"
"  -v                   verbose progress to stderr\n"
"  -h, --help           show this help\n"
"  --version            print version\n",
  DEFAULT_MAPQ, DEFAULT_MIN_LEN, DEFAULT_PRECISION, DEFAULT_THREADS);
}

static int parse_int(const char *flag, const char *s) {
  char *end = NULL;
  long v = strtol(s, &end, 10);
  if (end == NULL || *end != '\0') quit("Invalid integer for %s: %s", flag, s);
  return (int) v;
}

// HOMER makeTagDirectory flags that bam2td recognises only so it can warn
// and skip them — the user may have copied a makeTagDirectory command and
// pasted it as-is. n_args = number of value arguments to consume after the
// flag (-1 means "variadic", consume until the next - argument).
typedef struct homer_flag {
  const char *name;
  int         n_args;
} homer_flag_t;

static const homer_flag_t HOMER_IGNORED[] = {
  // bool flags
  { "-keep",                   0 }, { "-mask",                  0 },
  { "-update",                 0 }, { "-C",                     0 },
  { "-forceBED",               0 }, { "-assignMidPoint",        0 },
  { "-directional",            0 }, { "-force5th",              0 },
  { "-chrOnly",                0 }, { "-illuminaPE",            0 },
  { "-bowtiePE",               0 }, { "-checkGC",               0 },
  { "-removePEbg",             0 }, { "-removeSelfLigation",    0 },
  { "-removeRestrictionEnds",  0 }, { "-both",                  0 },
  { "-one",                    0 }, { "-onlyOne",               0 },
  { "-none",                   0 },
  // 1-arg flags
  { "-totalReads",             1 }, { "-normGC",                1 },
  { "-normLength",             1 }, { "-normOligo",             1 },
  { "-mCcontext",              1 }, { "-oligoStart",            1 },
  { "-oligoEnd",               1 }, { "-restrictionSiteLength", 1 },
  { "-minCounts",              1 }, { "-PEbgLength",            1 },
  { "-normFixedOligo",         1 }, { "-rsmis",                 1 },
  { "-freqStart",              1 }, { "-freqEnd",               1 },
  { "-iterNorm",               1 }, { "-minNormRatio",          1 },
  { "-maxNormRatio",           1 }, { "-restrictionSite",       1 },
  { "-name",                   1 }, { "-len",                   1 },
  { "-fragLength",             1 }, { "-format",                1 },
  // 2-arg flag
  { "-removeSpikes",           2 },
  // 3-arg flag
  { "-filterReads",            3 },
  // variadic
  { "-d",                     -1 }, { "-t",                    -1 },
};

#define HOMER_IGNORED_N (sizeof(HOMER_IGNORED) / sizeof(HOMER_IGNORED[0]))

// Return n_args (>=0), -1 for variadic, or -2 if `flag` is unknown.
static int homer_ignored_lookup(const char *flag) {
  for (size_t k = 0; k < HOMER_IGNORED_N; k++) {
    if (strcmp(flag, HOMER_IGNORED[k].name) == 0) return HOMER_IGNORED[k].n_args;
  }
  return -2;
}

static void parse_args(int argc, char **argv, opts_t *opts) {
  int i = 1;
  while (i < argc) {
    const char *a = argv[i];
    if (!strcmp(a, "-h") || !strcmp(a, "--help")) {
      usage(stdout); exit(EXIT_SUCCESS);
    } else if (!strcmp(a, "--version")) {
      puts("bam2td " BAM2TD_VERSION); exit(EXIT_SUCCESS);
    } else if (!strcmp(a, "-genome")) {
      if (++i >= argc) quit("-genome requires a name.");
      opts->genome = argv[i];
    } else if (!strcmp(a, "-mapq")) {
      if (++i >= argc) quit("-mapq requires an integer.");
      opts->mapq = parse_int(a, argv[i]);
    } else if (!strcmp(a, "-minlen")) {
      if (++i >= argc) quit("-minlen requires an integer.");
      opts->min_len = parse_int(a, argv[i]);
    } else if (!strcmp(a, "-maxlen")) {
      if (++i >= argc) quit("-maxlen requires an integer.");
      opts->max_len = parse_int(a, argv[i]);
    } else if (!strcmp(a, "-tbp")) {
      if (++i >= argc) quit("-tbp requires an integer.");
      opts->tbp = parse_int(a, argv[i]);
    } else if (!strcmp(a, "-precision")) {
      if (++i >= argc) quit("-precision requires 1, 2, or 3.");
      opts->precision = parse_int(a, argv[i]);
      if (opts->precision < 1 || opts->precision > 3) {
        quit("-precision must be 1, 2, or 3.");
      }
    } else if (!strcmp(a, "-mis")) {
      if (++i >= argc) quit("-mis requires an integer.");
      opts->mis = parse_int(a, argv[i]);
    } else if (!strcmp(a, "-p")) {
      if (++i >= argc) quit("-p requires an integer.");
      opts->threads = parse_int(a, argv[i]);
    } else if (!strcmp(a, "-flip"))       opts->flip       = true;
    else   if (!strcmp(a, "-sspe"))       opts->sspe       = true;
    else   if (!strcmp(a, "-single"))     opts->single     = true;
    else   if (!strcmp(a, "-rmsoft"))     opts->rmsoft     = true;
    else   if (!strcmp(a, "-omitSN"))     opts->omitSN     = true;
    else   if (!strcmp(a, "-keepAll"))    opts->keep_all   = true;
    else   if (!strcmp(a, "-keepOne"))    { opts->keep_one = true; opts->keep_all = false; }
    else   if (!strcmp(a, "-unique"))     opts->unique     = true;
    else   if (!strcmp(a, "-read1"))      opts->read1      = true;
    else   if (!strcmp(a, "-read2"))      opts->read2      = true;
    else   if (!strcmp(a, "-keepDup"))    opts->keep_dup   = true;
    else   if (!strcmp(a, "-keepQCfail")) opts->keep_qcfail = true;
    else   if (!strcmp(a, "-v"))          opts->v          = true;
    else if (a[0] == '-' && a[1] != '\0') {
      // Recognize HOMER makeTagDirectory flags we don't implement and skip
      // them with a warning, including any value arguments they consume.
      int n = homer_ignored_lookup(a);
      if (n >= 0) {
        warn("HOMER flag '%s' is not implemented by bam2td; ignoring.", a);
        for (int k = 0; k < n; k++) {
          if (i + 1 >= argc) break;
          i++;
        }
      } else if (n == -1) {
        warn("HOMER flag '%s' is not implemented by bam2td; ignoring its values too.", a);
        while (i + 1 < argc && argv[i + 1][0] != '-') i++;
      } else {
        quit("Unknown option: %s (try -h).", a);
      }
    } else {
      // Positionals, in order: <output-dir> <input.bam>.
      if (opts->output_dir == NULL)      opts->output_dir = a;
      else if (opts->input == NULL)      opts->input = a;
      else quit("Unexpected extra positional argument: '%s'.", a);
    }
    i++;
  }
  if (opts->output_dir == NULL)   quit("No output directory (try -h).");
  if (opts->input == NULL)        quit("No input file (try -h).");
  if (opts->read1 && opts->read2) quit("-read1 and -read2 are mutually exclusive.");
}

// Sorted-input check --------------------------------------------------------

static bool header_is_coord_sorted(bam_hdr_t *hdr) {
  const char *t = sam_hdr_str(hdr);
  if (t == NULL) return false;
  for (const char *p = t; *p != '\0'; ) {
    if (p[0] == '@' && p[1] == 'H' && p[2] == 'D' && (p[3] == '\t' || p[3] == ' ')) {
      const char *eol = strchr(p, '\n');
      if (eol == NULL) eol = p + strlen(p);
      for (const char *q = p; q + 3 <= eol; q++) {
        if (q[0] == 'S' && q[1] == 'O' && q[2] == ':') {
          const char *v = q + 3;
          const char *ve = v;
          while (ve < eol && *ve != '\t' && *ve != ' ' && *ve != '\n') ve++;
          return (ve - v == 10) && memcmp(v, "coordinate", 10) == 0;
        }
      }
      return false;
    }
    while (*p != '\0' && *p != '\n') p++;
    if (*p == '\n') p++;
  }
  return false;
}

// Per-read pipeline ---------------------------------------------------------

#define is_pri(b)         (((b)->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) == 0)
#define is_rev(b)         (((b)->core.flag & BAM_FREVERSE) != 0)
#define is_paired(b)      (((b)->core.flag & BAM_FPAIRED) != 0)
#define is_mate_unmap(b)  (((b)->core.flag & BAM_FMUNMAP) != 0)
#define is_proper_pair(b) (((b)->core.flag & BAM_FPROPER_PAIR) != 0)

typedef struct tag_t {
  int32_t   tid;
  hts_pos_t pos;       // 0-based: 5' end (SE) or fragment leftmost (PE)
  int8_t    strand;    // 0 = +, 1 = -
  int32_t   len;
  bool      is_pe_rep; // leftmost mate of a paired-end fragment
} tag_t;

static inline bool passes_filters(const bam1_t *b, const opts_t *opts) {
  const uint16_t f = b->core.flag;
  if (f & BAM_FUNMAP) return false;
  if ((f & BAM_FQCFAIL) && !opts->keep_qcfail) return false;
  if ((f & BAM_FDUP)    && !opts->keep_dup)    return false;
  if (!opts->keep_all && !is_pri(b))           return false;
  if (opts->read1 && !(f & BAM_FREAD1))        return false;
  if (opts->read2 && !(f & BAM_FREAD2))        return false;
  if (b->core.qual < opts->mapq)               return false;
  if (opts->unique) {
    const uint8_t *nh = bam_aux_get(b, "NH");
    if (nh != NULL && bam_aux2i(nh) != 1) return false;
  }
  if (opts->mis >= 0) {
    const uint8_t *nm = bam_aux_get(b, "NM");
    if (nm != NULL && bam_aux2i(nm) > opts->mis) return false;
  }
  return true;
}

// Read length excluding leading and trailing soft clips (mirrors quaqc's
// bases_covered() default branch).
static int32_t qlen_no_softclip(const bam1_t *b) {
  const int n_cigar = b->core.n_cigar;
  const uint32_t *cigar = bam_get_cigar(b);
  int32_t l = b->core.l_qseq;
  int kl, kr;
  for (kl = 0; kl < n_cigar; kl++) {
    if (bam_cigar_op(cigar[kl]) == BAM_CSOFT_CLIP) l -= bam_cigar_oplen(cigar[kl]);
    else break;
  }
  for (kr = n_cigar - 1; kr > kl; kr--) {
    if (bam_cigar_op(cigar[kr]) == BAM_CSOFT_CLIP) l -= bam_cigar_oplen(cigar[kr]);
    else break;
  }
  return l;
}

static bool compute_tag(const bam1_t *b, const opts_t *opts, tag_t *out) {
  const uint16_t f = b->core.flag;
  // A read is treated as part of a fragment only if it is a proper pair on
  // the same reference with a nonzero TLEN. Otherwise (chimeric, dovetailed,
  // mate-unmapped, mate-on-other-chr, |isize| absurdly large) fall back to
  // SE so we don't emit chromosome-spanning fragments.
  const bool pe_proper = is_paired(b) && is_proper_pair(b)
                         && !is_mate_unmap(b)
                         && (b->core.isize != 0)
                         && (b->core.tid == b->core.mtid);
  const bool pe_mode = pe_proper && !opts->single;

  // In PE mode each fragment must be counted once. Drop the rightmost mate
  // (isize < 0); only the leftmost mate (isize > 0) becomes a tag at b->pos.
  if (pe_mode && b->core.isize <= 0) return false;

  int32_t len;
  bool is_pe_rep = false;
  if (pe_mode) {
    len = (int32_t) b->core.isize;
    is_pe_rep = true;
  } else if (opts->rmsoft) {
    len = qlen_no_softclip(b);
  } else {
    len = b->core.l_qseq > 0 ? b->core.l_qseq
                              : (int32_t) bam_cigar2rlen(b->core.n_cigar, bam_get_cigar(b));
  }
  if (len < opts->min_len)                       return false;
  if (opts->max_len > 0 && len > opts->max_len)  return false;

  int strand = is_rev(b) ? 1 : 0;
  if (opts->sspe && pe_mode && (f & BAM_FREAD2)) strand ^= 1;
  if (opts->flip)                                strand ^= 1;

  // PE: fragment leftmost coord == b->core.pos (always sorted in coord-sorted
  // input). SE: 5' end of the read — for - strand reads that is bam_endpos-1,
  // which can land out of order vs adjacent + reads, so the chr stream sorts
  // per chromosome before emitting.
  hts_pos_t pos = b->core.pos;
  if (!pe_mode && is_rev(b)) {
    hts_pos_t end = bam_endpos(b);
    if (end > b->core.pos) pos = end - 1;
  }

  out->tid       = b->core.tid;
  out->pos       = pos;
  out->strand    = (int8_t) strand;
  out->len       = len;
  out->is_pe_rep = is_pe_rep;
  return true;
}

// Writers -------------------------------------------------------------------

typedef struct chr_stats_t {
  int64_t   unique_pos;
  double    total_tags;
  double    max_tbp;
  hts_pos_t max_pos_1based;   // HOMER's "appearentSize": max 1-based tag pos
                              // emitted to this chr's tag file. Sums to
                              // gsizeEstimate.
} chr_stats_t;

static char *sanitize_chr_filename(const char *name) {
  const size_t n = strlen(name);
  char *s = alloc(n + 1);
  bool warned = false;
  for (size_t i = 0; i < n; i++) {
    unsigned char c = (unsigned char) name[i];
    if (c == '/') {
      s[i] = '_';
      if (!warned) {
        warn("Replaced '/' in chr name '%s'.", name);
        warned = true;
      }
    } else {
      s[i] = (char) c;
    }
  }
  s[n] = '\0';
  return s;
}

// Histogram for medianTagsPerPosition. Counts are float but always
// integer-valued (each read = 1.0; -tbp caps but does not fractionate).
// Counts at or above MEDIAN_HIST_CAP fall into the overflow bucket; the
// median is still located correctly as long as it lies below the cap.
#define MEDIAN_HIST_CAP 4096
static int64_t median_hist[MEDIAN_HIST_CAP];

static inline void median_hist_record(double count) {
  int b = (int) (count + 0.5);
  if (b < 0) b = 0;
  if (b >= MEDIAN_HIST_CAP) b = MEDIAN_HIST_CAP - 1;
  median_hist[b]++;
}

static int median_hist_value(void) {
  int64_t total = 0;
  for (int i = 0; i < MEDIAN_HIST_CAP; i++) total += median_hist[i];
  if (total == 0) return 0;
  int64_t half = (total + 1) / 2;
  int64_t cum = 0;
  for (int i = 0; i < MEDIAN_HIST_CAP; i++) {
    cum += median_hist[i];
    if (cum >= half) return i;
  }
  return 0;
}

// chr_stream: buffers all tags for the current chromosome, sorts by
// (pos, strand, length) at the chromosome boundary, then emits one line per
// unique (pos, strand, length) tuple. Mirrors HOMER's makeTagDirectory,
// which never collapses across read lengths.
typedef struct chr_stream_t {
  const char  *out_dir;
  bam_hdr_t   *hdr;
  int          precision;
  int          tbp;            // 0 = unbounded
  chr_stats_t *stats;          // indexed by tid
  int          cur_tid;        // -1 = not started
  tag_t       *buf;
  size_t       n;
  size_t       cap;
} chr_stream_t;

static int tag_pos_strand_len_cmp(const void *a, const void *b) {
  const tag_t *x = (const tag_t *) a;
  const tag_t *y = (const tag_t *) b;
  if (x->pos    < y->pos)    return -1;
  if (x->pos    > y->pos)    return  1;
  if (x->strand != y->strand) return (int) x->strand - (int) y->strand;
  if (x->len    < y->len)    return -1;
  if (x->len    > y->len)    return  1;
  return 0;
}

static void chr_stream_flush(chr_stream_t *s) {
  if (s->n == 0 || s->cur_tid < 0) { s->n = 0; return; }
  qsort(s->buf, s->n, sizeof(tag_t), tag_pos_strand_len_cmp);

  const char *name = s->hdr->target_name[s->cur_tid];
  char *fname = sanitize_chr_filename(name);
  size_t need = strlen(s->out_dir) + 1 + strlen(fname) + sizeof(".tags.tsv");
  char *path = alloc(need);
  snprintf(path, need, "%s/%s.tags.tsv", s->out_dir, fname);
  FILE *fp = fopen(path, "w");
  if (fp == NULL) quit("Cannot open '%s' for writing: %s", path, strerror(errno));
  free(path);
  free(fname);

  chr_stats_t *cs = &s->stats[s->cur_tid];
  size_t i = 0;
  while (i < s->n) {
    const hts_pos_t pos    = s->buf[i].pos;
    const int8_t    strand = s->buf[i].strand;
    const int32_t   len    = s->buf[i].len;
    double          count  = 0.0;
    while (i < s->n
           && s->buf[i].pos    == pos
           && s->buf[i].strand == strand
           && s->buf[i].len    == len) {
      if (s->tbp <= 0 || count < (double) s->tbp) count += 1.0;
      i++;
    }
    fprintf(fp, "\t%s\t%" PRId64 "\t%d\t%.*f\t%d\n",
            name, (int64_t) (pos + 1), (int) strand,
            s->precision, count, (int) len);
    cs->unique_pos += 1;
    cs->total_tags += count;
    if (count > cs->max_tbp) cs->max_tbp = count;
    const hts_pos_t pos_1b = pos + 1;
    if (pos_1b > cs->max_pos_1based) cs->max_pos_1based = pos_1b;
    median_hist_record(count);
  }
  fclose(fp);
  s->n = 0;
}

static void chr_stream_emit(chr_stream_t *s, const tag_t *t) {
  if (s->cur_tid != t->tid) {
    if (s->cur_tid >= 0) chr_stream_flush(s);
    s->cur_tid = t->tid;
  }
  if (s->n == s->cap) {
    s->cap = s->cap ? s->cap * 2 : 65536;
    tag_t *nb = realloc(s->buf, s->cap * sizeof(tag_t));
    if (nb == NULL) quit("Out of memory (chr buffer %zu tags).", s->cap);
    s->buf = nb;
  }
  s->buf[s->n++] = *t;
}

static void chr_stream_close(chr_stream_t *s) {
  chr_stream_flush(s);
  free(s->buf);
  s->buf = NULL;
  s->cap = 0;
}

// pe_stream: petag.tsv is a single file across chromosomes, but we buffer per
// chromosome and sort by (pos, strand, length) before emitting so the line
// layout matches chr*.tags.tsv (one line per unique tuple, no length collapse).
typedef struct pe_stream_t {
  FILE       *fp;
  const char *out_dir;
  bam_hdr_t  *hdr;
  int         precision;
  int         tbp;
  int         cur_tid;
  tag_t      *buf;
  size_t      n;
  size_t      cap;
} pe_stream_t;

static void pe_stream_open_lazy(pe_stream_t *s) {
  if (s->fp != NULL) return;
  size_t need = strlen(s->out_dir) + sizeof("/petag.tsv");
  char *path = alloc(need);
  snprintf(path, need, "%s/petag.tsv", s->out_dir);
  s->fp = fopen(path, "w");
  if (s->fp == NULL) quit("Cannot open '%s' for writing: %s", path, strerror(errno));
  free(path);
}

static void pe_stream_flush(pe_stream_t *s) {
  if (s->n == 0 || s->cur_tid < 0) { s->n = 0; return; }
  qsort(s->buf, s->n, sizeof(tag_t), tag_pos_strand_len_cmp);
  pe_stream_open_lazy(s);

  const char *name = s->hdr->target_name[s->cur_tid];
  size_t i = 0;
  while (i < s->n) {
    const hts_pos_t pos    = s->buf[i].pos;
    const int8_t    strand = s->buf[i].strand;
    const int32_t   len    = s->buf[i].len;
    double          count  = 0.0;
    while (i < s->n
           && s->buf[i].pos    == pos
           && s->buf[i].strand == strand
           && s->buf[i].len    == len) {
      if (s->tbp <= 0 || count < (double) s->tbp) count += 1.0;
      i++;
    }
    fprintf(s->fp, "\t%s\t%" PRId64 "\t%d\t%.*f\t%d\n",
            name, (int64_t) (pos + 1), (int) strand,
            s->precision, count, (int) len);
  }
  s->n = 0;
}

static void pe_stream_emit(pe_stream_t *s, const tag_t *t) {
  if (s->cur_tid != t->tid) {
    if (s->cur_tid >= 0) pe_stream_flush(s);
    s->cur_tid = t->tid;
  }
  if (s->n == s->cap) {
    s->cap = s->cap ? s->cap * 2 : 65536;
    tag_t *nb = realloc(s->buf, s->cap * sizeof(tag_t));
    if (nb == NULL) quit("Out of memory (pe buffer %zu tags).", s->cap);
    s->buf = nb;
  }
  s->buf[s->n++] = *t;
}

static void pe_stream_close(pe_stream_t *s) {
  pe_stream_flush(s);
  free(s->buf);
  s->buf = NULL;
  s->cap = 0;
  if (s->fp != NULL) { fclose(s->fp); s->fp = NULL; }
}

// tagInfo.txt writer --------------------------------------------------------

static void write_tag_info(const opts_t *opts, bam_hdr_t *hdr,
                           chr_stats_t *stats, int n_targets,
                           bool saw_pe, double avg_tag_len,
                           int argc, char **argv)
{
  size_t need = strlen(opts->output_dir) + sizeof("/tagInfo.txt");
  char *path = alloc(need);
  snprintf(path, need, "%s/tagInfo.txt", opts->output_dir);
  FILE *fp = fopen(path, "w");
  if (fp == NULL) quit("Cannot open '%s': %s", path, strerror(errno));
  free(path);

  // gsizeEstimate mirrors HOMER: sum of the largest tag position emitted on
  // each chromosome (its `appearentSize`), not the full @SQ LN values from
  // the BAM header.
  int64_t total_unique = 0;
  double  total_tags   = 0.0;
  int64_t genome_size  = 0;
  for (int t = 0; t < n_targets; t++) {
    total_unique += stats[t].unique_pos;
    total_tags   += stats[t].total_tags;
    genome_size  += (int64_t) stats[t].max_pos_1based;
  }
  const double tags_per_bp = genome_size > 0
                               ? total_tags / (double) genome_size : 0.0;
  const double avg_per_pos = total_unique > 0
                               ? total_tags / (double) total_unique : 0.0;
  const int    median_tpp  = median_hist_value();

  // HOMER's layout: row 1 is a 3-column header, row 2 is the summary row
  // (genome name, total unique, total tags), then key=value metadata, then
  // 3-column per-chromosome rows, then a cmd= trailer. HOMER prefixes
  // `genome=` to the value only when `-genome <name>` was supplied;
  // otherwise the first cell is the bare literal `genome`.
  fprintf(fp, "name\tUnique Positions\tTotal Tags\n");
  if (opts->genome != NULL) {
    fprintf(fp, "genome=%s\t%" PRId64 "\t%.1f\n",
            opts->genome, total_unique, total_tags);
  } else {
    fprintf(fp, "genome\t%" PRId64 "\t%.1f\n", total_unique, total_tags);
  }

  // HOMER pads every single-cell line with two trailing tabs so the file
  // visually keeps three "columns". We match that padding to stay diff-clean
  // against makeTagDirectory output.
  fprintf(fp, "fragmentLengthEstimate=%d\t\t\n",
          saw_pe ? 0 : opts->default_fraglen);
  fprintf(fp, "peakSizeEstimate=%d\t\t\n",
          saw_pe ? 0 : opts->default_fraglen);
  fprintf(fp, "tagsPerBP=%.6f\t\t\n", tags_per_bp);
  fprintf(fp, "averageTagsPerPosition=%.3f\t\t\n", avg_per_pos);
  fprintf(fp, "medianTagsPerPosition=%d\t\t\n", median_tpp);
  fprintf(fp, "averageTagLength=%.3f\t\t\n", avg_tag_len);
  fprintf(fp, "gsizeEstimate=%" PRId64 "\t\t\n", genome_size);
  fprintf(fp, "averageFragmentGCcontent=%.3f\t\t\n", saw_pe ? 0.0 : -1.0);

  for (int t = 0; t < n_targets; t++) {
    if (opts->omitSN && stats[t].unique_pos == 0) continue;
    fprintf(fp, "%s\t%" PRId64 "\t%.1f\n",
            hdr->target_name[t],
            stats[t].unique_pos,
            stats[t].total_tags);
  }

  fputs("cmd=", fp);
  for (int i = 0; i < argc; i++) {
    if (i > 0) fputc(' ', fp);
    fputs(argv[i], fp);
  }
  fputs("\t\t\n", fp);
  fclose(fp);
}

// Main loop -----------------------------------------------------------------

int main(int argc, char **argv) {
  const time_t time_start = time(NULL);
  setlocale(LC_NUMERIC, "en_US.UTF-8");   // enable %' thousands separators

  opts_t opts_storage = default_opts();
  opts_t *opts = &opts_storage;
  parse_args(argc, argv, opts);

  if (mkdir(opts->output_dir, 0777) < 0 && errno != EEXIST) {
    quit("Cannot create directory '%s': %s", opts->output_dir, strerror(errno));
  }

  htsFile *bam = hts_open(opts->input, "r");
  if (bam == NULL) quit("Failed to open '%s'.", opts->input);
  if (opts->threads > 1) hts_set_threads(bam, opts->threads);

  bam_hdr_t *hdr = sam_hdr_read(bam);
  if (hdr == NULL) quit("Failed to read header from '%s'.", opts->input);
  if (!header_is_coord_sorted(hdr)) {
    quit("Input must be coordinate-sorted (@HD SO:coordinate).");
  }
  if (hdr->n_targets <= 0) quit("No reference sequences in header.");

  const int n_targets = hdr->n_targets;
  chr_stats_t *cstats = alloc(sizeof(chr_stats_t) * (size_t) n_targets);

  chr_stream_t chr_stream = {0};
  chr_stream.out_dir   = opts->output_dir;
  chr_stream.hdr       = hdr;
  chr_stream.precision = opts->precision;
  chr_stream.tbp       = opts->tbp;
  chr_stream.stats     = cstats;
  chr_stream.cur_tid   = -1;

  pe_stream_t pe_stream = {0};
  pe_stream.out_dir   = opts->output_dir;
  pe_stream.hdr       = hdr;
  pe_stream.precision = opts->precision;
  pe_stream.tbp       = opts->tbp;

  bam1_t *aln = bam_init1();
  if (aln == NULL) quit("Out of memory.");

  bool    saw_pe         = false;
  int64_t global_len_sum = 0;
  int64_t global_len_n   = 0;
  int     ret;

  msg("Reading %s\n", opts->input);

  while ((ret = sam_read1(bam, hdr, aln)) >= 0) {
    if (!passes_filters(aln, opts)) continue;
    tag_t tag;
    if (!compute_tag(aln, opts, &tag)) continue;

    chr_stream_emit(&chr_stream, &tag);
    if (tag.is_pe_rep && !opts->single) {
      pe_stream_emit(&pe_stream, &tag);
      saw_pe = true;
    }
    global_len_sum += tag.len;
    global_len_n   += 1;
  }
  if (ret < -1) warn("Truncated input.");

  chr_stream_close(&chr_stream);
  pe_stream_close(&pe_stream);

  const double avg_len = global_len_n > 0
                           ? (double) global_len_sum / (double) global_len_n
                           : 0.0;
  write_tag_info(opts, hdr, cstats, n_targets, saw_pe, avg_len, argc, argv);

  msg("Wrote tag directory '%s' (%'" PRId64 " tags%s).\n",
      opts->output_dir, global_len_n, saw_pe ? ", paired-end" : "");

  free(cstats);
  bam_destroy1(aln);
  sam_hdr_destroy(hdr);
  hts_close(bam);

  if (opts->v) {
    fprintf(stderr, "Time elapsed: ");
    print_time(stderr, (time_t) difftime(time(NULL), time_start));
    print_mem(stderr);
  }
  return EXIT_SUCCESS;
}
