# bam2td

Standalone replacement for HOMER's `makeTagDirectory` when the input is a
BAM/SAM/CRAM file. Reads coordinate-sorted alignments via htslib and writes
a tag directory that downstream HOMER tools (`findPeaks`, `annotatePeaks.pl`,
`analyzeRepeats.pl`, ...) can consume directly.

The implementation is single-pass and never holds the whole file in RAM.
Tags are buffered per chromosome, sorted by `(pos, strand, length)` at each
chromosome boundary, and emitted as collapsed runs — one line per unique
`(pos, strand, length)` tuple, matching HOMER's `makeTagDirectory` output
layout exactly.

It is **not** intended to be a feature-complete replacement; only the bare
essential files are created. HOMER `tagAutocorrelation.txt`,
`tagLengthDistribution.txt`, `tagCountDistribution.txt`, and `tagGCcontent.txt`
files are not produced. Run `makeTagDirectory <tagdir> -udpate -checkGC -genome <genome>`
afterwards to complete the tag directory.

## Build

Vendored htslib and zlib source trees live under `libs/`, are built as static
archives, and linked into the final binary.

```sh
# Default release build (statically links libs/htslib/libhts.a and libs/zlib/libz.a)
make release

# Build the bundled libraries from source first, then link
make release-full

# Debug build with -fsanitize=address,undefined
make debug

# Link against the system htslib / zlib instead of bundled sources
make release hts_dyn=1 z_dyn=1

# Native CPU targeting
make release native=1

# Optional htslib features (off by default)
make release-full with_curl=1 with_bz2=1 with_lzma=1
```

`make release` expects `libs/htslib/libhts.a` and `libs/zlib/libz.a` to
already exist. Use `make release-full` to build them from the bundled
sources, or pass `hts_dyn=1` / `z_dyn=1` to use the system libraries.

## Usage

```sh
bam2td <tagdir> <input.bam|sam|cram> [options]
```

The input must be coordinate-sorted (`@HD SO:coordinate` in the header);
otherwise the tool refuses to run. Program arguments mirror that of
`makeTagDirectory`.

### Options

```
-genome <name>       genome name written into tagInfo.txt (default unknown)
-mapq <int>          minimum MAPQ (default 10)
-minlen <int>        minimum tag/fragment length (default 0)
-maxlen <int>        maximum tag/fragment length (0 = no cap)
-tbp <int>           cap reads per (pos, strand, length) (0 = no cap)
-precision <1|2|3>   decimal places when printing counts (default 1)
-mis <int>           maximum NM:i (mismatches+indels); -1 disables
-flip                flip strand of all reads
-sspe                strand-specific PE: flip strand of R2
-single              treat input as single-end (skip petag.tsv)
-rmsoft              report query length excluding soft-clipped bases
-omitSN              omit unused @SQ entries from tagInfo.txt
-keepAll             keep secondary and supplementary alignments
-keepOne             keep primary alignments only (the default)
-unique              require NH:i==1 where present
-read1 / -read2      in PE data, keep only R1 or R2
-keepDup             do not skip BAM_FDUP reads (default: skip)
-keepQCfail          do not skip BAM_FQCFAIL reads (default: skip)
-p <int>             htslib decoder threads (default 1)
-v                   verbose progress + final elapsed time and peak memory
-h, --help           show help
--version            print version
```

### Behaviour notes

- **Paired-end handling.** When the alignment is flagged
  `BAM_FPROPER_PAIR | BAM_FPAIRED` and the mate maps to the same reference
  with a nonzero TLEN, the leftmost mate is emitted as a single fragment
  tag (`pos = b->core.pos`, `len = abs(isize)`, `strand = strand of the
  leftmost mate`), and the rightmost mate is dropped. Improper pairs,
  chimeric reads, and mate-unmapped reads fall back to single-end treatment
  (each emitted at its own 5' end). `-single` forces SE treatment for the
  entire file.
- **Filters.** By default `BAM_FUNMAP`, `BAM_FQCFAIL`, `BAM_FDUP`, and
  non-primary alignments (`BAM_FSECONDARY | BAM_FSUPPLEMENTARY`) are
  skipped. `-keepDup`, `-keepQCfail`, `-keepAll` opt back in. `-mapq`
  filters on `b->core.qual`; `-mis` reads `NM:i`; `-unique` requires
  `NH:i==1` when present.
- **Length.** For SE reads the default tag length is the read length
  (`b->core.l_qseq`). `-rmsoft` drops leading and trailing soft-clipped
  bases. For PE-rep tags the length is `abs(isize)` and is unaffected by
  `-rmsoft`.
- **Strand.** Forward reads are strand `0`, reverse `1`. `-flip` inverts
  every read; `-sspe` additionally inverts R2 inside a proper pair.
- **`-tbp`** caps the number of reads collapsed into one
  `(pos, strand, length)` bucket. The cap applies to both `chr*.tags.tsv`
  and `petag.tsv`.

## Output

```
<tagdir>/<chr>.tags.tsv     5' tag records, sorted by position
<tagdir>/petag.tsv          fragment records (PE only)
<tagdir>/tagInfo.txt        summary header consumed by HOMER
```

Each `<chr>.tags.tsv` and `petag.tsv` line has a leading tab (HOMER's
convention so the chromosome column lands in field 2) followed by five
fields:

```
\t<chr>\t<1-based pos>\t<strand>\t<count>\t<length>
```

One line per unique `(pos, strand, length)` tuple, sorted by position
within each chromosome.

`tagInfo.txt` mirrors HOMER's layout:

1. Header row: `name\tUnique Positions\tTotal Tags`
2. Summary row: `<genome>\t<total unique>\t<total tags>` (the first cell is
   `genome=<name>` if `-genome` was supplied, otherwise the bare literal
   `genome`).
3. `key=value` metadata lines for `fragmentLengthEstimate`, `peakSizeEstimate`,
   `tagsPerBP`, `averageTagsPerPosition`, `medianTagsPerPosition`,
   `averageTagLength`, `gsizeEstimate`, `averageFragmentGCcontent`. Each is
   right-padded with two trailing tabs to keep three "columns".
4. Per-chromosome rows: `<chr>\t<unique positions>\t<total tags>`.
5. `cmd=<full command line>` trailer.

`gsizeEstimate` mirrors HOMER's `appearentSize` accounting: the sum of the
largest tag position emitted on each chromosome (i.e. the right edge of the
data), not the sum of `@SQ LN:` lengths.

## Limitations

- Requires coordinate-sorted input. Pipe through `samtools sort` if
  needed.
- CRAM input needs a reference; configure `REF_PATH`/`REF_CACHE` per
  htslib conventions.
- Memory cost is one `tag_t` (~32 bytes, including buffer-doubling
  overhead) per pre-collapse passing read on the *largest* chromosome —
  the buffer is freed and reallocated on each chromosome boundary, so the
  whole genome never sits in RAM. As a rough yardstick, ~63M tags on the
  deepest chromosome corresponds to ~2 GB peak RSS; typical ChIP-seq /
  ATAC-seq / csRNA-seq workloads land well under 100 MB.
