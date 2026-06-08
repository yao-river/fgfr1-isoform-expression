#!/bin/bash
# ============================================================
# Extract FGFR1 transcript TPM from all salmon quant.sf on GCS
#
# RUN IN: Cloud Shell, or a VM in us-central1 (SAME region as the bucket,
#         so reads incur no egress charge). Only the small TSV is taken out.
#
# Cost model: streams each quant.sf via `gsutil cat` and greps FGFR1 rows
# only — full files are never saved to disk. Same-region GCS->compute reads
# are free; the result TSV is a few MB.
# ============================================================
set -uo pipefail

# ---- CONFIG ----
BUCKET="gs://salmon-quant"           # per-sample folders live here
PROJECT="sellers-lab-yaohe-dev"      # billing project for -u
ENST_FILE="fgfr1_ensts.txt"          # one ENST id per line (upload alongside)
OUT="fgfr1_transcript_tpm.tsv"
QUANT_NAME="quant.sf"
PARALLEL=32                          # concurrent gsutil cat workers
# ----------------

# Build grep pattern: bare ENST id followed by '.'(version) or TAB,
# so it matches whether or not quant.sf keeps the .version suffix.
PATTERN=$(paste -sd'|' "$ENST_FILE")
GREP_RE="^(${PATTERN})[.$(printf '\t')]"
export PROJECT QUANT_NAME GREP_RE

# Worker: stream one sample's quant.sf, keep only FGFR1 rows, prefix sample id
extract_one() {
    local sample_path="$1"
    local name; name=$(basename "${sample_path%/}")
    gsutil -u "$PROJECT" cat "${sample_path}${QUANT_NAME}" 2>/dev/null \
        | grep -E "$GREP_RE" \
        | awk -v s="$name" 'BEGIN{OFS="\t"} {print s, $0}'
}
export -f extract_one

# Header
echo -e "sample\tName\tLength\tEffectiveLength\tTPM\tNumReads" > "$OUT"

# List sample folders, fan out in parallel, append results
gsutil -u "$PROJECT" ls "${BUCKET}/" \
    | xargs -P "$PARALLEL" -I{} bash -c 'extract_one "$@"' _ {} \
    >> "$OUT"

echo "Done. Rows: $(($(wc -l < "$OUT") - 1)); output: $OUT" >&2