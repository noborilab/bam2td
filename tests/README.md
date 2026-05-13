# tests

Synthetic test suite for `bam2td`. All test BAMs are constructed on the fly
from inline SAM text and piped through `samtools sort`, so this directory
stays tiny and clean to commit.

## Running

```sh
make release        # build the binary first
make test           # runs bash tests/run.sh
```

Or directly:

```sh
bash tests/run.sh
```

To exercise a custom binary (e.g. a debug build):

```sh
BAM2TD=./bam2td-debug bash tests/run.sh
```

## Dependencies

- `samtools` on `PATH`
- `bash`, `awk`, `grep`, `sed`, `xxd` (standard on macOS / any Linux)

## Coverage

The suite verifies argument handling, error paths, output format
(leading-tab line shape, per-`(pos, strand, length)` line splitting, sort
order, 1-based positions), paired-end semantics (proper-pair filtering,
isize-derived length, `-single` fallback, `-read1` selection), filter
flags (`-mapq`, `-minlen`, `-maxlen`, `-tbp`, `-mis`, `-unique`,
`-keepDup`, `-keepAll`), strand and length transforms (`-flip`,
`-rmsoft`), and the HOMER-compatible `tagInfo.txt` layout (column header
row, summary row with optional `genome=` prefix, three-column per-chr
rows, `appearentSize`-based `gsizeEstimate`, `-omitSN`, `-precision`).
