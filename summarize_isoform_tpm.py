"""
Summarize FGFR1 transcript-level TPM into isoform-level TPM.

Inputs:
  - fgfr1_transcript_tpm.tsv     : output of extract_fgfr1_tpm.sh
        cols: sample, Name, Length, EffectiveLength, TPM, NumReads
  - fgfr1_isoform_classified.csv : transcript_id -> isoform mapping

Output:
  - fgfr1_isoform_tpm.csv        : isoform-level TPM matrix (sample x isoform)

Notes:
  - Only samples with the full set of transcripts (EXPECTED_N) are kept;
    samples with fewer rows (incomplete/truncated reads) are dropped.
  - salmon TPM is already length-corrected and globally normalized over the
    full transcriptome. We only SUM transcript TPMs within each isoform class.
    We do NOT renormalize or re-correct length.
"""
import pandas as pd

TX_TPM      = "./results/fgfr1_transcript_tpm.tsv"
ISOFORM_MAP = "./results/fgfr1_isoform_classified.csv"
OUT         = "./results/fgfr1_isoform_tpm.csv"
EXPECTED_N  = 73     # keep only samples with the full set of FGFR1 transcripts

# ---- load transcript-level TPM (long format) ----
tx = pd.read_csv(TX_TPM, sep="\t")

# ---- drop incomplete samples ----
counts = tx.groupby("sample").size()
full_samples = counts[counts == EXPECTED_N].index
n_dropped = tx["sample"].nunique() - len(full_samples)
tx = tx[tx["sample"].isin(full_samples)].copy()
print(f"Kept {len(full_samples)} complete samples "
      f"(== {EXPECTED_N} transcripts); dropped {n_dropped} incomplete.")

# ---- load isoform map ----
iso = pd.read_csv(ISOFORM_MAP)[["transcript_id", "isoform"]]

# strip version suffix on both sides so IDs match regardless of salmon index
tx["enst"]  = tx["Name"].str.split(".").str[0]
iso["enst"] = iso["transcript_id"].str.split(".").str[0]

# ---- join transcript -> isoform ----
merged = tx.merge(iso[["enst", "isoform"]], on="enst", how="left")
n_unmapped = merged["isoform"].isna().sum()
if n_unmapped:
    print(f"[warn] {n_unmapped} rows had no isoform mapping:")
    print(merged[merged["isoform"].isna()]["Name"].unique()[:10])

# ---- sum TPM within (sample, isoform), pivot to sample x isoform ----
iso_tpm = (merged.dropna(subset=["isoform"])
           .groupby(["sample", "isoform"])["TPM"].sum()
           .unstack(fill_value=0.0))

iso_tpm.to_csv(OUT)
print(f"\nOutput: {OUT}")
print(f"Matrix: {iso_tpm.shape[0]} samples x {iso_tpm.shape[1]} isoforms")
print("\nIsoform columns:")
print(iso_tpm.columns.tolist())

# ---- quick sanity: mean TPM per isoform across samples ----
print("\nMean TPM per isoform (across samples):")
print(iso_tpm.mean().sort_values(ascending=False).round(2).to_string())