# In [1]: imports and paths
from pathlib import Path

import pandas as pd


BASE = Path(__file__).resolve().parent
HUMAN_PATH = BASE / "circRNA_presence_master.v4.tsv"
MOUSE_PATH = BASE / "mouse_circRNA_presence_master.v4.hg38_1based.tsv"
OUT_PATH = BASE / "circRNA_presence_master.v4.mouse_conservation.tsv"


# In [2]: helper for strand normalization
NULL_STRINGS = {"", "nan", "na", "none", "null", "<na>"}


def normalize_strand_series(series: pd.Series) -> pd.Series:
    text = series.astype("string").str.strip()
    lowered = text.str.lower()
    return text.mask(text.isna() | lowered.isin(NULL_STRINGS), "NA")


# In [3]: read human hg38_1based circRNA presence table
human_df = pd.read_csv(
    HUMAN_PATH,
    sep="\t",
    usecols=["chrom", "start", "end", "strand"],
    dtype={"chrom": str, "start": "int64", "end": "int64", "strand": "string"},
    keep_default_na=True,
)


# In [4]: normalize human strand and build human coordinate ID
human_df["strand"] = normalize_strand_series(human_df["strand"])
human_df["circ_id"] = (
    human_df["chrom"]
    + ":"
    + human_df["start"].astype(str)
    + "|"
    + human_df["end"].astype(str)
    + "("
    + human_df["strand"].astype(str)
    + ")"
)


# In [5]: read mouse circRNAs lifted from mm10_1based to hg38_1based
mouse_df = pd.read_csv(
    MOUSE_PATH,
    sep="\t",
    usecols=[
        "chrom",
        "start",
        "end",
        "strand",
        "original_chrom",
        "original_start",
        "original_end",
        "original_strand",
    ],
    dtype={
        "chrom": str,
        "start": "int64",
        "end": "int64",
        "strand": "string",
        "original_chrom": str,
        "original_start": "int64",
        "original_end": "int64",
        "original_strand": "string",
    },
    keep_default_na=True,
)


# In [6]: normalize mouse strand columns and build coordinate IDs
mouse_df["strand"] = normalize_strand_series(mouse_df["strand"])
mouse_df["original_strand"] = normalize_strand_series(mouse_df["original_strand"])

mouse_df["circ_id"] = (
    mouse_df["chrom"]
    + ":"
    + mouse_df["start"].astype(str)
    + "|"
    + mouse_df["end"].astype(str)
    + "("
    + mouse_df["strand"].astype(str)
    + ")"
)

mouse_df["mouse_circ_id"] = (
    mouse_df["original_chrom"]
    + ":"
    + mouse_df["original_start"].astype(str)
    + "|"
    + mouse_df["original_end"].astype(str)
    + "("
    + mouse_df["original_strand"].astype(str)
    + ")"
)


# In [7]: collapse multiple mouse circRNAs mapped to the same human coordinate
mouse_df_grouped = (
    mouse_df[["circ_id", "mouse_circ_id"]]
    .drop_duplicates()
    .groupby("circ_id", sort=False)["mouse_circ_id"]
    .agg(lambda values: ",".join(values))
    .reset_index()
)


# In [8]: left join mouse conservation evidence onto human circRNAs
merged_df = human_df.merge(mouse_df_grouped, on="circ_id", how="left")
merged_df["mouse_circ_id"] = merged_df["mouse_circ_id"].fillna("")


# In [9]: write output
merged_df.to_csv(OUT_PATH, sep="\t", index=False)


# In [10]: validation summary
print(f"human_rows={len(human_df)}")
print(f"human_NA_strand_rows={(human_df['strand'] == 'NA').sum()}")
print(f"mouse_lifted_rows={len(mouse_df)}")
print(f"mouse_NA_strand_rows={(mouse_df['strand'] == 'NA').sum()}")
print(f"mouse_NA_original_strand_rows={(mouse_df['original_strand'] == 'NA').sum()}")
print(f"mouse_grouped_circ_ids={len(mouse_df_grouped)}")
print(f"output_rows={len(merged_df)}")
print(f"conserved_human_rows={(merged_df['mouse_circ_id'] != '').sum()}")
print(
    "conserved_NA_strand_rows="
    f"{((merged_df['strand'] == 'NA') & (merged_df['mouse_circ_id'] != '')).sum()}"
)
