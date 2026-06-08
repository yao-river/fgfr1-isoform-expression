"""
FGFR1 isoform classification from GENCODE annotation.

Input : GENCODE GTF subset for FGFR1 (with transcript/exon structure)
Output: per-transcript anchor coverage matrix + isoform class label

Design:
  - Driven by anchor exon coordinates; switching genes only requires editing CONFIG.
  - Uses fractional "coverage" instead of binary within-overlap, so truncated
    exons (partial coverage of an anchor) are handled correctly.
  - Classification decision tree is based on structural completeness
    (D3 domain / mutually exclusive IIIb-IIIc exons / kinase domain).
"""
import pandas as pd
import re

# ============================================================
# CONFIG -- edit only this block when switching to another gene
# ============================================================
GENE_NAME = "FGFR1"
CHROM = "chr8"

ANCHORS = {  # (start, end), 1-based, GRCh38
    "IgI":       (38429682, 38429948),   # alpha/beta switch (exon 3 / IgI)
    "IgIII_pre": (38424509, 38424685),   # invariant first half of D3 (IIIa region)
    "IIIb":      (38423035, 38423175),   # mutually exclusive IIIb exon
    "IIIc":      (38421798, 38421941),   # full IIIc exon (baseline for truncation)
    "TK_1": (38413918, 38414022), "TK_2": (38414151, 38414288),
    "TK_3": (38414558, 38414629), "TK_4": (38414779, 38414901),
    "TK_5": (38415870, 38416061), "TK_6": (38417307, 38417414),
    "TK_7": (38417874, 38417990),
}
TK_KEYS = [k for k in ANCHORS if k.startswith("TK_")]

# Thresholds (set after inspecting the real coverage distribution)
TH = {
    "anchor_hit":     0.5,   # a single anchor with coverage >= this counts as "hit"
    "iiic_full":      0.9,   # IIIc >= this -> full
    "iiic_trunc_min": 0.3,   # IIIc in [trunc_min, full) -> truncated
    "iiib_full":      0.9,
    "tk_full":        7,     # number of TK segments hit >= this -> full kinase domain
    "d3_entered":     0.5,   # IgIII_pre >= this -> transcript enters D3
    "igI_present":    0.5,   # IgI >= this -> alpha
}

# ============================================================
# 1. Parse GTF (pure pandas, no third-party dependency)
# ============================================================
def parse_gtf(path):
    cols = ["seqname", "source", "feature", "start", "end",
            "score", "strand", "frame", "attr"]
    df = pd.read_csv(path, sep="\t", header=None, names=cols, comment="#")

    def attr(a, k):
        m = re.search(rf'{k} "([^"]+)"', a)
        return m.group(1) if m else None

    df["transcript_id"]   = df["attr"].apply(lambda a: attr(a, "transcript_id"))
    df["transcript_type"] = df["attr"].apply(lambda a: attr(a, "transcript_type"))
    return df

# ============================================================
# 2. Anchor coverage matrix
# ============================================================
def overlap_bp(es, ee, a_s, a_e):
    return max(0, min(ee, a_e) - max(es, a_s) + 1)

def build_matrix(df):
    exons = df[df["feature"] == "exon"]
    tx_bio = (df[df["feature"] == "transcript"]
              .set_index("transcript_id")["transcript_type"].to_dict())
    rows = []
    for tx, g in exons.groupby("transcript_id"):
        rec = {"transcript_id": tx, "biotype": tx_bio.get(tx, "NA")}
        ex = list(zip(g["start"], g["end"]))
        for name, (a_s, a_e) in ANCHORS.items():
            alen = a_e - a_s + 1
            # total overlap of all exons with this anchor, as a fraction of anchor length
            tot = sum(overlap_bp(es, ee, a_s, a_e) for es, ee in ex)
            rec[name] = round(tot / alen, 3)
        rec["n_exons"] = len(ex)
        rows.append(rec)
    mat = pd.DataFrame(rows).set_index("transcript_id")
    # number of TK segments covered above the hit threshold
    mat["TK_hits"] = (mat[TK_KEYS] >= TH["anchor_hit"]).sum(axis=1)
    return mat

# ============================================================
# 3. Classification decision tree
# ============================================================
def classify(r):
    if r["biotype"] != "protein_coding":
        return "Noncoding"

    has_igI = r["IgI"]       >= TH["igI_present"]
    in_d3   = r["IgIII_pre"] >= TH["d3_entered"]
    iiib    = r["IIIb"]      >= TH["iiib_full"]
    iiic_full  = r["IIIc"] >= TH["iiic_full"]
    iiic_trunc = TH["iiic_trunc_min"] <= r["IIIc"] < TH["iiic_full"]
    tk_full = r["TK_hits"]   >= TH["tk_full"]

    ab = "alpha" if has_igI else "beta"

    # ---- Full kinase domain: full-length signaling receptor ----
    if tk_full:
        if iiic_full:   return f"{ab}_IIIc"
        if iiic_trunc:  return f"{ab}_IIIc_truncated"
        if iiib:        return f"{ab}_IIIb"
        return "coding_other"

    # ---- Absent / partial kinase domain ----
    # Still has IIIb/IIIc but incomplete kinase domain -> dominant-negative candidate
    if iiib or iiic_full or iiic_trunc:
        return "dominant_negative_candidate"
    # Enters D3 but no IIIb/IIIc and no kinase domain -> secreted IIIa candidate
    if in_d3:
        return f"{ab}_IIIa_secreted_candidate"
    # Truncated before reaching D3:
    #   has IgI            -> alpha_truncated
    #   no IgI (all other anchors also 0, since D3 not reached) -> ultrashort
    if has_igI:
        return "alpha_truncated"
    return "ultrashort_truncated"

# ============================================================
# main
# ============================================================
def run(gtf_path, out_path):
    df = parse_gtf(gtf_path)
    mat = build_matrix(df)
    mat["isoform"] = mat.apply(classify, axis=1)
    order = (["biotype", "isoform", "IgI", "IgIII_pre", "IIIb", "IIIc"]
             + TK_KEYS + ["TK_hits", "n_exons"])
    mat = mat[order].sort_values(["isoform"])
    mat.to_csv(out_path)
    return mat

if __name__ == "__main__":
    import sys
    gtf = sys.argv[1] if len(sys.argv) > 1 else "/tmp/fgfr1_v49.gtf"
    out = sys.argv[2] if len(sys.argv) > 2 else "/tmp/fgfr1_isoform_classified.csv"
    mat = run(gtf, out)
    print(f"=== {GENE_NAME} isoform classification ({len(mat)} transcripts) ===\n")
    print(mat["isoform"].value_counts().to_string())
    print("\n=== protein_coding detail ===")
    pd.set_option("display.width", 200, "display.max_columns", 30)
    pc = mat[mat["biotype"] == "protein_coding"]
    print(pc[["isoform", "IgI", "IgIII_pre", "IIIb", "IIIc", "TK_hits"]].to_string())
    print(f"\nOutput: {out}")