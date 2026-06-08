#!/bin/bash
# ============================================================
# Extract FGFR1 transcript TPM from all salmon quant.sf on GCS
#
# RUN IN: Cloud Shell, or a VM in us-central1 (SAME region as the bucket,
#         so reads incur no egress charge).
#
# FIX vs previous version: each worker writes its OWN temp file, then we
# concatenate at the end. This avoids the parallel-append interleaving that
# corrupted some lines (xargs -P writing to one file is NOT atomic above 4KB).
# ============================================================
set -uo pipefail

# ---- CONFIG ----
BUCKET="gs://salmon-quant"
PROJECT="sellers-lab-yaohe-dev"
ENST_FILE="./results/fgfr1_ensts.txt"          # one bare ENST id per line
OUT="./results/fgfr1_transcript_tpm.tsv"
QUANT_NAME="quant.sf"
PARALLEL=32
TMPDIR_X="fgfr1_tpm_tmp"             # per-sample temp files go here
# ----------------

# grep pattern: bare ENST followed by '.'(version) or TAB, matches anywhere
# in the line (ENST may be in col 1 OR col 2 depending on the salmon index).
PATTERN=$(paste -sd'|' "$ENST_FILE")
GREP_RE="(${PATTERN})[.$(printf '\t')]"
export PROJECT QUANT_NAME GREP_RE TMPDIR_X

mkdir -p "$TMPDIR_X"

# Worker: stream one sample, keep FGFR1 rows, write to its OWN file.
extract_one() {
    local sample_path="$1"
    local name; name=$(basename "${sample_path%/}")
    local out="${TMPDIR_X}/${name}.tsv"
    # awk prints sample id + the original line, guarantees trailing newline
    gsutil -u "$PROJECT" cat "${sample_path}${QUANT_NAME}" 2>/dev/null \
        | grep -E "$GREP_RE" \
        | awk -v s="$name" 'BEGIN{OFS="\t"} {print s, $0}' \
        > "$out"
    # remove empty result (sample with no quant.sf / read error)
    [ -s "$out" ] || rm -f "$out"
}
export -f extract_one

echo "Listing samples..." >&2
gsutil -u "$PROJECT" ls "${BUCKET}/" \
    | xargs -P "$PARALLEL" -I{} bash -c 'extract_one "$@"' _ {}

echo "Merging per-sample files..." >&2
# header + concatenate all per-sample files
{
    echo -e "sample\tName\tLength\tEffectiveLength\tTPM\tNumReads"
    cat "${TMPDIR_X}"/*.tsv
} > "$OUT"

n_samples=$(ls "${TMPDIR_X}" | wc -l)
n_rows=$(($(wc -l < "$OUT") - 1))
echo "Done. Samples: $n_samples ; rows: $n_rows ; output: $OUT" >&2
echo "Temp files kept in ${TMPDIR_X}/ (rm -rf it after you verify)." >&2