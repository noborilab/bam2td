# bam2td

`bam2td` is a small standalone tool that stands in for HOMER's
`makeTagDirectory` when the input you happen to have is a BAM/SAM/CRAM file.
It reads coordinate-sorted alignments through htslib and writes out a tag
directory that the downstream HOMER tools (`findPeaks`, `annotatePeaks.pl`,
`analyzeRepeats.pl`, and so on) can consume directly.

It works in a single pass and never holds the whole file in memory at once.
Tags are buffered one chromosome at a time, sorted by `(pos, strand, length)`
at each chromosome boundary, and then written out as collapsed runs (one line
per unique `(pos, strand, length)` tuple), which is exactly the layout
`makeTagDirectory` produces.

It isn't meant to be a feature-complete replacement, to be honest; only the
bare essential files get created. In particular the HOMER
`tagAutocorrelation.txt`, `tagLengthDistribution.txt`,
`tagCountDistribution.txt`, and `tagGCcontent.txt` files are not produced. If
you do need those, running `makeTagDirectory <tagdir> -update -checkGC -genome <genome>`
afterwards will fill the rest of the tag directory in.

## Build

The htslib and zlib source trees are vendored under `libs/`; they get built as
static archives and linked straight into the final binary.

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

`make release` expects `libs/htslib/libhts.a` and `libs/zlib/libz.a` to be
there already. If they aren't, `make release-full` builds them from the
bundled sources first, and passing `hts_dyn=1` / `z_dyn=1` links against the
system libraries instead.

## Usage

```sh
bam2td <tagdir> <input.bam|sam|cram> [options]
```

The input has to be coordinate-sorted (`@HD SO:coordinate` in the header); if
it isn't, the tool refuses to run rather than quietly produce something wrong.
The argument order and the option names mirror `makeTagDirectory`, so for the
most part an existing command should carry over.

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

- **Paired-end handling.** When an alignment is flagged
  `BAM_FPROPER_PAIR | BAM_FPAIRED` and its mate maps to the same reference
  with a nonzero TLEN, the leftmost mate is emitted as a single fragment tag
  (`pos = b->core.pos`, `len = abs(isize)`, with the strand taken from the
  leftmost mate) and the rightmost mate is dropped. Improper pairs, chimeric
  reads, and mate-unmapped reads all fall back to single-end treatment (each
  emitted at its own 5' end). `-single` forces SE treatment for the whole
  file.
- **Filters.** By default `BAM_FUNMAP`, `BAM_FQCFAIL`, `BAM_FDUP`, and the
  non-primary alignments (`BAM_FSECONDARY | BAM_FSUPPLEMENTARY`) are skipped;
  `-keepDup`, `-keepQCfail`, and `-keepAll` opt back in. `-mapq` filters on
  `b->core.qual`, `-mis` reads `NM:i`, and `-unique` requires `NH:i==1`
  wherever it is present.
- **Length.** For SE reads the default tag length is the read length
  (`b->core.l_qseq`), and `-rmsoft` drops the leading and trailing
  soft-clipped bases. For PE-rep tags the length is `abs(isize)`, which
  `-rmsoft` does not touch.
- **Strand.** Forward reads are strand `0` and reverse reads strand `1`.
  `-flip` inverts every read, and `-sspe` additionally inverts R2 inside a
  proper pair.
- **`-tbp`** caps how many reads get collapsed into a single
  `(pos, strand, length)` bucket (the cap applies to both `chr*.tags.tsv` and
  `petag.tsv`).

## Output

```
<tagdir>/<chr>.tags.tsv     5' tag records, sorted by position
<tagdir>/petag.tsv          fragment records (PE only)
<tagdir>/tagInfo.txt        summary header consumed by HOMER
```

Each `<chr>.tags.tsv` and `petag.tsv` line starts with a leading tab (this is
HOMER's convention, so that the chromosome column lands in field 2), followed
by five fields:

```
\t<chr>\t<1-based pos>\t<strand>\t<count>\t<length>
```

There is one line per unique `(pos, strand, length)` tuple, sorted by position
within each chromosome.

`tagInfo.txt` mirrors HOMER's layout:

1. Header row: `name\tUnique Positions\tTotal Tags`.
2. Summary row: `<genome>\t<total unique>\t<total tags>` (the first cell is
   `genome=<name>` if `-genome` was supplied, otherwise just the bare literal
   `genome`).
3. `key=value` metadata lines for `fragmentLengthEstimate`,
   `peakSizeEstimate`, `tagsPerBP`, `averageTagsPerPosition`,
   `medianTagsPerPosition`, `averageTagLength`, `gsizeEstimate`, and
   `averageFragmentGCcontent`. Each one is right-padded with two trailing tabs
   to keep the three "columns" lined up.
4. Per-chromosome rows: `<chr>\t<unique positions>\t<total tags>`.
5. A `cmd=<full command line>` trailer.

One thing worth noting (since it is easy to assume otherwise):
`gsizeEstimate` follows HOMER's `appearentSize` accounting, which is the sum
of the largest tag position emitted on each chromosome (the right edge of the
data, in other words), rather than the sum of the `@SQ LN:` lengths.

## Limitations

- The input has to be coordinate-sorted; pipe it through `samtools sort`
  first if it isn't.
- CRAM input needs a reference, so you will want to set `REF_PATH` /
  `REF_CACHE` per the usual htslib conventions.
- Memory cost is roughly one `tag_t` (~32 bytes, including the buffer-doubling
  overhead) per pre-collapse passing read on the *largest* chromosome. Because
  the buffer is freed and reallocated at each chromosome boundary, the whole
  genome never sits in RAM at once. As a rough yardstick, ~63M tags on the
  deepest chromosome works out to ~2 GB peak RSS, while typical ChIP-seq,
  ATAC-seq, and csRNA-seq workloads land comfortably under 100 MB.
