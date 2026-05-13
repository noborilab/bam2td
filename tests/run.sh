#!/bin/bash
# bam2td synthetic test suite.
#
# Requires `samtools` on PATH (used to build tiny coordinate-sorted BAMs
# from inline SAM text). The bam2td binary is taken from $BAM2TD if set,
# otherwise the binary one directory above this script.
#
# Run from anywhere:
#   bash tests/run.sh
# Or:
#   make test

set -uo pipefail

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
BAM2TD=${BAM2TD:-"$(cd "$SCRIPT_DIR/.." && pwd)/bam2td"}

if [ ! -x "$BAM2TD" ]; then
  echo "bam2td binary not found at $BAM2TD — run 'make release' first." >&2
  exit 2
fi
if ! command -v samtools >/dev/null 2>&1; then
  echo "samtools not on PATH — required to build test BAMs." >&2
  exit 2
fi

TMPROOT=$(mktemp -d /tmp/bam2td_tests_XXXXX)
PASS=0; FAIL=0
trap 'rm -rf "$TMPROOT"' EXIT

# ---------------------------------------------------------------------------
# Test framework
# ---------------------------------------------------------------------------

run_test() {
  local name="$1"
  local d="$TMPROOT/$name"
  mkdir -p "$d"
  printf 'Testing %-45s ' "$name ..."
  if "test_$name" "$d" >"$d/stdout" 2>"$d/stderr"; then
    echo "PASS"
    PASS=$((PASS + 1))
  else
    echo "FAIL"
    sed 's/^/    /' "$d/stderr" | tail -5 >&2
    FAIL=$((FAIL + 1))
  fi
}

assert_eq() {
  local label="$1" actual="$2" expected="$3"
  if [ "$actual" = "$expected" ]; then return 0; fi
  echo "  assert_eq $label: expected='$expected' got='$actual'" >&2
  return 1
}

assert_file_exists() {
  if [ -s "$1" ]; then return 0; fi
  echo "  assert_file_exists: '$1' missing or empty" >&2
  return 1
}

assert_file_absent() {
  if [ ! -e "$1" ]; then return 0; fi
  echo "  assert_file_absent: '$1' unexpectedly present" >&2
  return 1
}

assert_contains() {
  local file="$1" pattern="$2"
  if grep -qF "$pattern" "$file"; then return 0; fi
  echo "  assert_contains: '$pattern' not in $file" >&2
  return 1
}

assert_not_contains() {
  local file="$1" pattern="$2"
  if ! grep -qF "$pattern" "$file"; then return 0; fi
  echo "  assert_not_contains: '$pattern' unexpectedly in $file" >&2
  return 1
}

assert_lines() {
  local file="$1" expected="$2"
  local n; n=$(wc -l <"$file" | tr -d ' ')
  assert_eq "lines($file)" "$n" "$expected"
}

# ---------------------------------------------------------------------------
# SAM emission helpers (positions are 1-based, SAM convention)
# ---------------------------------------------------------------------------

# Repeat a single character N times to fill SEQ.
_seq() { printf '%*s' "$1" '' | tr ' ' "$2"; }

# chr_header chr1 LEN [chr2 LEN ...]
chr_header() {
  printf '@HD\tVN:1.6\tSO:queryname\n'
  while [ $# -ge 2 ]; do
    printf '@SQ\tSN:%s\tLN:%d\n' "$1" "$2"
    shift 2
  done
}

# se_read name chr pos len mapq [flag] [extra_aux]
se_read() {
  local name="$1" chr="$2" pos="$3" len="$4" mapq="${5:-42}" flag="${6:-0}"
  local aux="${7:-}"
  local s; s=$(_seq "$len" A)
  printf '%s\t%d\t%s\t%d\t%d\t%dM\t*\t0\t0\t%s\t*' \
    "$name" "$flag" "$chr" "$pos" "$mapq" "$len" "$s"
  if [ -n "$aux" ]; then printf '\t%s' "$aux"; fi
  printf '\n'
}

# pe_pair name chr r1_pos r1_len r2_pos r2_len [mapq] [extra_r1_flag] [extra_r2_flag]
# R1 is +strand, R2 is -strand, proper pair.
pe_pair() {
  local name="$1" chr="$2" r1p="$3" r1l="$4" r2p="$5" r2l="$6"
  local mapq="${7:-42}" xf1="${8:-0}" xf2="${9:-0}"
  local tlen=$((r2p + r2l - r1p))
  local f1=$(( 99 + xf1 ))    # paired+proper+mate_rev+r1
  local f2=$(( 147 + xf2 ))   # paired+proper+rev+r2
  local s1 s2
  s1=$(_seq "$r1l" A); s2=$(_seq "$r2l" T)
  printf '%s\t%d\t%s\t%d\t%d\t%dM\t=\t%d\t%d\t%s\t*\n' \
    "$name" "$f1" "$chr" "$r1p" "$mapq" "$r1l" "$r2p" "$tlen" "$s1"
  printf '%s\t%d\t%s\t%d\t%d\t%dM\t=\t%d\t-%d\t%s\t*\n' \
    "$name" "$f2" "$chr" "$r2p" "$mapq" "$r2l" "$r1p" "$tlen" "$s2"
}

# Compile SAM on stdin to a coordinate-sorted, indexed BAM.
make_bam() {
  samtools sort -o "$1" - 2>/dev/null
}

# ---------------------------------------------------------------------------
# Tests — argument handling
# ---------------------------------------------------------------------------

test_help_exits_zero() {
  "$BAM2TD" -h >/dev/null
}

test_version_exits_zero() {
  local out; out=$("$BAM2TD" --version)
  echo "$out" | grep -q '^bam2td '
}

test_missing_args_errors() {
  if "$BAM2TD" 2>/dev/null; then return 1; fi
  return 0
}

test_old_o_flag_rejected() {
  if "$BAM2TD" -o /tmp/x foo.bam 2>/dev/null; then return 1; fi
  return 0
}

test_unknown_flag_rejected() {
  local d="$1"
  { chr_header chr1 1000; se_read r1 chr1 100 50; } | make_bam "$d/in.bam"
  if "$BAM2TD" --not-a-real-flag "$d/td" "$d/in.bam" 2>/dev/null; then return 1; fi
  return 0
}

test_read1_read2_mutex() {
  local d="$1"
  { chr_header chr1 1000; se_read r1 chr1 100 50; } | make_bam "$d/in.bam"
  if "$BAM2TD" -read1 -read2 "$d/td" "$d/in.bam" 2>/dev/null; then return 1; fi
  return 0
}

test_unsorted_input_rejected() {
  local d="$1"
  # Write a BAM that the header marks as queryname-sorted, not coordinate.
  { chr_header chr1 1000; se_read r1 chr1 100 50; } \
    | samtools view -bS - 2>/dev/null > "$d/in.bam"
  if "$BAM2TD" "$d/td" "$d/in.bam" 2>/dev/null; then
    echo "  expected non-zero exit on unsorted input" >&2
    return 1
  fi
  return 0
}

test_missing_input_errors() {
  local d="$1"
  if "$BAM2TD" "$d/td" "$d/does_not_exist.bam" 2>/dev/null; then return 1; fi
  return 0
}

# ---------------------------------------------------------------------------
# Tests — basic SE output
# ---------------------------------------------------------------------------

test_se_basic_outputs() {
  local d="$1"
  {
    chr_header chr1 1000
    se_read r1 chr1 100 50
    se_read r2 chr1 200 50
    se_read r3 chr1 300 50
  } | make_bam "$d/in.bam"
  "$BAM2TD" "$d/td" "$d/in.bam" >/dev/null
  assert_file_exists "$d/td/chr1.tags.tsv" || return 1
  assert_file_exists "$d/td/tagInfo.txt"   || return 1
  assert_file_absent "$d/td/petag.tsv"     || return 1
  assert_lines       "$d/td/chr1.tags.tsv" 3
}

test_se_positions_1based_in_output() {
  local d="$1"
  { chr_header chr1 1000; se_read r1 chr1 100 50; } | make_bam "$d/in.bam"
  "$BAM2TD" "$d/td" "$d/in.bam" >/dev/null
  # SAM pos is 1-based -> bam2td emits 1-based -> column 2 must be 100.
  local p; p=$(awk '{print $2}' "$d/td/chr1.tags.tsv")
  assert_eq "tag pos" "$p" "100"
}

test_se_leading_tab() {
  local d="$1"
  { chr_header chr1 1000; se_read r1 chr1 100 50; } | make_bam "$d/in.bam"
  "$BAM2TD" "$d/td" "$d/in.bam" >/dev/null
  # First byte of the file must be a tab (HOMER's convention).
  local first; first=$(head -c 1 "$d/td/chr1.tags.tsv" | xxd -p)
  assert_eq "first byte (hex)" "$first" "09"
}

test_se_reverse_strand_5prime() {
  local d="$1"
  # A 100bp + read at pos 200 ends at 299; a 100bp - read with leftmost
  # pos 50 also ends at 149; 5' pos for - is 149 (1-based).
  {
    chr_header chr1 1000
    se_read fwd chr1 200 100 42 0
    se_read rev chr1 50  100 42 16
  } | make_bam "$d/in.bam"
  "$BAM2TD" "$d/td" "$d/in.bam" >/dev/null
  # Both lines should appear; - strand 5' is 149.
  assert_contains "$d/td/chr1.tags.tsv" $'\t200\t0\t1.0\t100' || return 1
  assert_contains "$d/td/chr1.tags.tsv" $'\t149\t1\t1.0\t100'
}

test_se_lines_sorted_by_pos() {
  local d="$1"
  {
    chr_header chr1 10000
    # Mix +/- in non-monotonic 5' order; output must be sorted by pos.
    se_read a chr1 100 50 42 0     # 5' = 100
    se_read b chr1  50 50 42 16    # 5' = 99
    se_read c chr1 300 50 42 0     # 5' = 300
    se_read d chr1 250 50 42 16    # 5' = 299
  } | make_bam "$d/in.bam"
  "$BAM2TD" "$d/td" "$d/in.bam" >/dev/null
  local sorted; sorted=$(sort -k2,2n "$d/td/chr1.tags.tsv")
  local actual; actual=$(cat "$d/td/chr1.tags.tsv")
  [ "$sorted" = "$actual" ] || { echo "  output not sorted by pos"; return 1; }
}

test_se_collapse_by_pos_strand_length() {
  local d="$1"
  # 5 reads all at the same 5' (pos 100, + strand) but two distinct lengths
  # must collapse into TWO lines, not one (HOMER-style length split).
  {
    chr_header chr1 1000
    for i in 1 2 3; do se_read "r${i}" chr1 100 50; done
    for i in 4 5;   do se_read "r${i}" chr1 100 75; done
  } | make_bam "$d/in.bam"
  "$BAM2TD" "$d/td" "$d/in.bam" >/dev/null
  assert_lines "$d/td/chr1.tags.tsv" 2 || return 1
  # Verify the count column sums correctly per length.
  local c50 c75
  c50=$(awk '$2==100 && $3==0 && $5==50 {print $4}' "$d/td/chr1.tags.tsv")
  c75=$(awk '$2==100 && $3==0 && $5==75 {print $4}' "$d/td/chr1.tags.tsv")
  assert_eq "count len=50"  "$c50" "3.0" || return 1
  assert_eq "count len=75"  "$c75" "2.0"
}

# ---------------------------------------------------------------------------
# Tests — PE handling
# ---------------------------------------------------------------------------

test_pe_basic_outputs() {
  local d="$1"
  {
    chr_header chr1 5000
    pe_pair p1 chr1 100 50 300 50
    pe_pair p2 chr1 500 50 700 50
  } | make_bam "$d/in.bam"
  "$BAM2TD" "$d/td" "$d/in.bam" >/dev/null
  assert_file_exists "$d/td/chr1.tags.tsv" || return 1
  assert_file_exists "$d/td/petag.tsv"     || return 1
  # One tag per fragment (leftmost mate only) -> 2 lines, not 4.
  assert_lines "$d/td/chr1.tags.tsv" 2
}

test_pe_length_is_isize() {
  local d="$1"
  # R1 at 100, R2 at 300 with 50bp reads: fragment spans 100..349 -> 250bp.
  {
    chr_header chr1 5000
    pe_pair p1 chr1 100 50 300 50
  } | make_bam "$d/in.bam"
  "$BAM2TD" "$d/td" "$d/in.bam" >/dev/null
  local len; len=$(awk '{print $5}' "$d/td/chr1.tags.tsv")
  assert_eq "PE tag len = isize" "$len" "250"
}

test_pe_single_flag_disables_petag() {
  local d="$1"
  {
    chr_header chr1 5000
    pe_pair p1 chr1 100 50 300 50
    pe_pair p2 chr1 500 50 700 50
  } | make_bam "$d/in.bam"
  "$BAM2TD" -single "$d/td" "$d/in.bam" >/dev/null
  assert_file_absent "$d/td/petag.tsv" || return 1
  # In SE mode both mates emit; expect 4 lines.
  assert_lines "$d/td/chr1.tags.tsv" 4
}

test_pe_improper_falls_back_to_se() {
  local d="$1"
  # Same as a PE pair but clear the proper-pair bit (0x2). FLAG 99 -> 97,
  # FLAG 147 -> 145.
  {
    chr_header chr1 5000
    printf 'p\t97\tchr1\t100\t42\t50M\t=\t300\t250\t%s\t*\n' "$(_seq 50 A)"
    printf 'p\t145\tchr1\t300\t42\t50M\t=\t100\t-250\t%s\t*\n' "$(_seq 50 T)"
  } | make_bam "$d/in.bam"
  "$BAM2TD" "$d/td" "$d/in.bam" >/dev/null
  # Both mates should now be emitted as independent SE tags.
  assert_lines "$d/td/chr1.tags.tsv" 2 || return 1
  assert_file_absent "$d/td/petag.tsv"
}

test_pe_read1_filter() {
  local d="$1"
  {
    chr_header chr1 5000
    pe_pair p1 chr1 100 50 300 50
  } | make_bam "$d/in.bam"
  # -read1 + -single: only R1 (the +-strand mate at pos 100).
  "$BAM2TD" -single -read1 "$d/td" "$d/in.bam" >/dev/null
  assert_lines "$d/td/chr1.tags.tsv" 1 || return 1
  assert_contains "$d/td/chr1.tags.tsv" $'\t100\t0\t'
}

# ---------------------------------------------------------------------------
# Tests — filter flags
# ---------------------------------------------------------------------------

test_mapq_filter() {
  local d="$1"
  {
    chr_header chr1 1000
    se_read hi chr1 100 50 42
    se_read lo chr1 200 50 5
  } | make_bam "$d/in.bam"
  "$BAM2TD" -mapq 30 "$d/td" "$d/in.bam" >/dev/null
  assert_lines "$d/td/chr1.tags.tsv" 1
}

test_minlen_filter() {
  local d="$1"
  {
    chr_header chr1 1000
    se_read short chr1 100 25
    se_read long  chr1 200 75
  } | make_bam "$d/in.bam"
  "$BAM2TD" -minlen 50 "$d/td" "$d/in.bam" >/dev/null
  assert_lines "$d/td/chr1.tags.tsv" 1 || return 1
  assert_contains "$d/td/chr1.tags.tsv" $'\t200\t'
}

test_maxlen_filter() {
  local d="$1"
  {
    chr_header chr1 1000
    se_read short chr1 100 25
    se_read long  chr1 200 75
  } | make_bam "$d/in.bam"
  "$BAM2TD" -maxlen 50 "$d/td" "$d/in.bam" >/dev/null
  assert_lines "$d/td/chr1.tags.tsv" 1 || return 1
  assert_contains "$d/td/chr1.tags.tsv" $'\t100\t'
}

test_tbp_caps_count() {
  local d="$1"
  {
    chr_header chr1 1000
    for i in 1 2 3 4 5; do se_read "r${i}" chr1 100 50; done
  } | make_bam "$d/in.bam"
  "$BAM2TD" -tbp 1 "$d/td" "$d/in.bam" >/dev/null
  local c; c=$(awk '{print $4}' "$d/td/chr1.tags.tsv")
  assert_eq "tbp=1 count" "$c" "1.0"
}

test_mis_filter() {
  local d="$1"
  {
    chr_header chr1 1000
    se_read ok  chr1 100 50 42 0 'NM:i:0'
    se_read bad chr1 200 50 42 0 'NM:i:5'
  } | make_bam "$d/in.bam"
  "$BAM2TD" -mis 2 "$d/td" "$d/in.bam" >/dev/null
  assert_lines "$d/td/chr1.tags.tsv" 1 || return 1
  assert_contains "$d/td/chr1.tags.tsv" $'\t100\t'
}

test_unique_filter() {
  local d="$1"
  {
    chr_header chr1 1000
    se_read uniq chr1 100 50 42 0 'NH:i:1'
    se_read mult chr1 200 50 42 0 'NH:i:3'
  } | make_bam "$d/in.bam"
  "$BAM2TD" -unique "$d/td" "$d/in.bam" >/dev/null
  assert_lines "$d/td/chr1.tags.tsv" 1 || return 1
  assert_contains "$d/td/chr1.tags.tsv" $'\t100\t'
}

test_dup_skipped_by_default() {
  local d="$1"
  {
    chr_header chr1 1000
    se_read ok  chr1 100 50 42 0       # FLAG 0
    se_read dup chr1 200 50 42 1024    # FLAG 0x400 = DUP
  } | make_bam "$d/in.bam"
  "$BAM2TD" "$d/td" "$d/in.bam" >/dev/null
  assert_lines "$d/td/chr1.tags.tsv" 1
}

test_keepDup_keeps_duplicates() {
  local d="$1"
  {
    chr_header chr1 1000
    se_read ok  chr1 100 50 42 0
    se_read dup chr1 200 50 42 1024
  } | make_bam "$d/in.bam"
  "$BAM2TD" -keepDup "$d/td" "$d/in.bam" >/dev/null
  assert_lines "$d/td/chr1.tags.tsv" 2
}

test_secondary_skipped_by_default() {
  local d="$1"
  {
    chr_header chr1 1000
    se_read pri chr1 100 50 42 0       # primary
    se_read sec chr1 200 50 42 256     # 0x100 secondary
  } | make_bam "$d/in.bam"
  "$BAM2TD" "$d/td" "$d/in.bam" >/dev/null
  assert_lines "$d/td/chr1.tags.tsv" 1
}

test_keepAll_keeps_secondary() {
  local d="$1"
  {
    chr_header chr1 1000
    se_read pri chr1 100 50 42 0
    se_read sec chr1 200 50 42 256
  } | make_bam "$d/in.bam"
  "$BAM2TD" -keepAll "$d/td" "$d/in.bam" >/dev/null
  assert_lines "$d/td/chr1.tags.tsv" 2
}

# ---------------------------------------------------------------------------
# Tests — strand and length flags
# ---------------------------------------------------------------------------

test_flip_inverts_strand() {
  local d="$1"
  { chr_header chr1 1000; se_read r chr1 100 50 42 0; } | make_bam "$d/in.bam"
  "$BAM2TD" -flip "$d/td" "$d/in.bam" >/dev/null
  # FWD read with -flip should appear as strand=1.
  local s; s=$(awk '{print $3}' "$d/td/chr1.tags.tsv")
  assert_eq "flipped strand" "$s" "1"
}

test_rmsoft_excludes_softclip() {
  local d="$1"
  # CIGAR 10S40M -> soft-clipped read, 40bp on reference, l_qseq=50.
  # Without -rmsoft -> length 50; with -rmsoft -> length 40.
  {
    chr_header chr1 1000
    printf 'r\t0\tchr1\t100\t42\t10S40M\t*\t0\t0\t%s\t*\n' "$(_seq 50 A)"
  } | make_bam "$d/in.bam"
  "$BAM2TD" "$d/td" "$d/in.bam" >/dev/null
  local lwithout; lwithout=$(awk '{print $5}' "$d/td/chr1.tags.tsv")
  assert_eq "no -rmsoft length" "$lwithout" "50" || return 1

  rm -rf "$d/td2"
  "$BAM2TD" -rmsoft "$d/td2" "$d/in.bam" >/dev/null
  local lwith; lwith=$(awk '{print $5}' "$d/td2/chr1.tags.tsv")
  assert_eq "-rmsoft length" "$lwith" "40"
}

# ---------------------------------------------------------------------------
# Tests — tagInfo.txt format
# ---------------------------------------------------------------------------

test_taginfo_header_layout() {
  local d="$1"
  {
    chr_header chr1 1000
    se_read r chr1 100 50
  } | make_bam "$d/in.bam"
  "$BAM2TD" "$d/td" "$d/in.bam" >/dev/null
  # Line 1 column headers (HOMER convention).
  local l1; l1=$(sed -n '1p' "$d/td/tagInfo.txt")
  assert_eq "header line" "$l1" $'name\tUnique Positions\tTotal Tags' || return 1
  # Line 2 is summary row with bare 'genome' literal (no -genome).
  local l2; l2=$(sed -n '2p' "$d/td/tagInfo.txt")
  echo "$l2" | grep -qE '^genome\t[0-9]+\t[0-9.]+$' \
    || { echo "  summary row: '$l2'"; return 1; }
  # Has required metadata keys.
  assert_contains "$d/td/tagInfo.txt" 'fragmentLengthEstimate=' || return 1
  assert_contains "$d/td/tagInfo.txt" 'medianTagsPerPosition='  || return 1
  assert_contains "$d/td/tagInfo.txt" 'gsizeEstimate='          || return 1
  assert_contains "$d/td/tagInfo.txt" 'cmd='
}

test_taginfo_genome_flag() {
  local d="$1"
  { chr_header chr1 1000; se_read r chr1 100 50; } | make_bam "$d/in.bam"
  "$BAM2TD" -genome hg99 "$d/td" "$d/in.bam" >/dev/null
  # Summary row starts with `genome=hg99` (only when -genome supplied).
  assert_contains "$d/td/tagInfo.txt" 'genome=hg99'
}

test_taginfo_perchr_three_columns() {
  local d="$1"
  {
    chr_header chr1 1000 chr2 1000
    se_read a chr1 100 50
    se_read b chr2 200 50
  } | make_bam "$d/in.bam"
  "$BAM2TD" "$d/td" "$d/in.bam" >/dev/null
  # Per-chr rows must be 3 tab-separated columns: chr\tunique\ttotal.
  local row; row=$(grep -E '^chr1\b' "$d/td/tagInfo.txt")
  local cols; cols=$(awk -F'\t' '{print NF}' <<<"$row")
  assert_eq "per-chr columns" "$cols" "3" || return 1
  echo "$row" | grep -qE '^chr1\t[0-9]+\t[0-9.]+$' \
    || { echo "  bad row '$row'"; return 1; }
}

test_taginfo_gsize_uses_max_pos() {
  local d="$1"
  # chr1 has length 10000 in header but the rightmost tag is at pos 200
  # (1-based). gsizeEstimate should be 200, NOT 10000 (HOMER's appearentSize).
  {
    chr_header chr1 10000
    se_read a chr1 100 50
    se_read b chr1 200 50
  } | make_bam "$d/in.bam"
  "$BAM2TD" "$d/td" "$d/in.bam" >/dev/null
  local gs; gs=$(grep -oE 'gsizeEstimate=[0-9]+' "$d/td/tagInfo.txt" \
                  | cut -d= -f2)
  assert_eq "gsizeEstimate" "$gs" "200"
}

test_omitSN_skips_empty_chrs() {
  local d="$1"
  {
    chr_header chr1 1000 chr2 1000
    se_read a chr1 100 50          # chr2 has zero tags
  } | make_bam "$d/in.bam"
  "$BAM2TD" -omitSN "$d/td" "$d/in.bam" >/dev/null
  # chr2 should not appear in tagInfo.txt when -omitSN is set.
  if grep -qE '^chr2\b' "$d/td/tagInfo.txt"; then
    echo "  chr2 unexpectedly present in tagInfo.txt"
    return 1
  fi
}

test_precision_decimals() {
  local d="$1"
  { chr_header chr1 1000; se_read r chr1 100 50; } | make_bam "$d/in.bam"
  "$BAM2TD" -precision 3 "$d/td" "$d/in.bam" >/dev/null
  # Count column should have three decimals: "1.000".
  local c; c=$(awk '{print $4}' "$d/td/chr1.tags.tsv")
  assert_eq "precision=3 count" "$c" "1.000"
}

# ---------------------------------------------------------------------------
# Tests — multi-chromosome sanity
# ---------------------------------------------------------------------------

test_multi_chr_separate_files() {
  local d="$1"
  {
    chr_header chr1 1000 chr2 1000
    se_read a chr1 100 50
    se_read b chr2 200 50
  } | make_bam "$d/in.bam"
  "$BAM2TD" "$d/td" "$d/in.bam" >/dev/null
  assert_file_exists "$d/td/chr1.tags.tsv" || return 1
  assert_file_exists "$d/td/chr2.tags.tsv" || return 1
  assert_lines "$d/td/chr1.tags.tsv" 1     || return 1
  assert_lines "$d/td/chr2.tags.tsv" 1
}

# ---------------------------------------------------------------------------
# Run all tests
# ---------------------------------------------------------------------------

ALL=(
  # argument handling
  help_exits_zero
  version_exits_zero
  missing_args_errors
  old_o_flag_rejected
  unknown_flag_rejected
  read1_read2_mutex
  unsorted_input_rejected
  missing_input_errors
  # SE basics
  se_basic_outputs
  se_positions_1based_in_output
  se_leading_tab
  se_reverse_strand_5prime
  se_lines_sorted_by_pos
  se_collapse_by_pos_strand_length
  # PE basics
  pe_basic_outputs
  pe_length_is_isize
  pe_single_flag_disables_petag
  pe_improper_falls_back_to_se
  pe_read1_filter
  # filters
  mapq_filter
  minlen_filter
  maxlen_filter
  tbp_caps_count
  mis_filter
  unique_filter
  dup_skipped_by_default
  keepDup_keeps_duplicates
  secondary_skipped_by_default
  keepAll_keeps_secondary
  # strand / length
  flip_inverts_strand
  rmsoft_excludes_softclip
  # tagInfo.txt
  taginfo_header_layout
  taginfo_genome_flag
  taginfo_perchr_three_columns
  taginfo_gsize_uses_max_pos
  omitSN_skips_empty_chrs
  precision_decimals
  # multi-chr
  multi_chr_separate_files
)

for t in "${ALL[@]}"; do run_test "$t"; done

echo
echo "================================"
echo "Results: $PASS passed, $FAIL failed out of $((PASS + FAIL)) tests"
echo

[ "$FAIL" -eq 0 ]
