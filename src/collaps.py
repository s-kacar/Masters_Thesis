
import time, pandas as pd
from pathlib import Path          # ← or: import pathlib



# Load your melted matrix
df = pd.read_csv("E:\Guido\sibel\Masters_Thesis\data\hog_intensity_long.tsv", sep='\t')
df = df.copy()  # avoid SettingWithCopyWarning
# Identify all LFQ columns
lfq_cols = [col for col in df.columns if 'MaxLFQ' in col]  # adjust as needed
print(lfq_cols)
# Load your melted matrix
df = pd.read_csv("E:\Guido\sibel\Masters_Thesis\data\hog_intensity_long.tsv", sep='\t')
df = df.copy()

# List of tissue columns (adjust these to match your actual column names)
lfq_cols = [
    'P095141_1 MaxLFQ Intensity', 'P095142_2 MaxLFQ Intensity', 'P095143_3 MaxLFQ Intensity',
    'P095144_4 MaxLFQ Intensity', 'P095145_5 MaxLFQ Intensity', 'P095146_6 MaxLFQ Intensity',
    'P095147_7 MaxLFQ Intensity', 'P095148_8 MaxLFQ Intensity', 'P095149_9 MaxLFQ Intensity',
    'P095150_10 MaxLFQ Intensity', 'P095151_11 MaxLFQ Intensity', 'P095152_12 MaxLFQ Intensity',
    'P095153_13 MaxLFQ Intensity', 'P095154_14 MaxLFQ Intensity', 'P095155_15 MaxLFQ Intensity',
    'P0104005_1 MaxLFQ Intensity', 'P090741_2 MaxLFQ Intensity', 'P090742_3 MaxLFQ Intensity',
    'P090743_4 MaxLFQ Intensity', 'P090744_5 MaxLFQ Intensity', 'P090745_6 MaxLFQ Intensity',
    'P090746_7 MaxLFQ Intensity', 'P090747_8 MaxLFQ Intensity', 'P090748_9 MaxLFQ Intensity',
    'P090749_10 MaxLFQ Intensity', 'P090750_11 MaxLFQ Intensity', 'P090751_12 MaxLFQ Intensity',
    'P090752_13 MaxLFQ Intensity', 'P090753_14 MaxLFQ Intensity', 'P090754_15 MaxLFQ Intensity',
    'P090755_16 MaxLFQ Intensity', 'P090756_17 MaxLFQ Intensity', 'P092517_1 MaxLFQ Intensity',
    'P092518_2 MaxLFQ Intensity', 'P092519_3 MaxLFQ Intensity', 'P092520_4 MaxLFQ Intensity',
    'P092521_5 MaxLFQ Intensity', 'P092522_6 MaxLFQ Intensity', 'P092523_7 MaxLFQ Intensity',
    'P092524_8 MaxLFQ Intensity', 'P092525_9 MaxLFQ Intensity', 'P092526_10 MaxLFQ Intensity',
    'P092527_11 MaxLFQ Intensity', 'P092528_12 MaxLFQ Intensity', 'P092529_13 MaxLFQ Intensity',
    'P092530_14 MaxLFQ Intensity', 'P092531_15 MaxLFQ Intensity', 'P092532_16 MaxLFQ Intensity',
    'P091307_1 MaxLFQ Intensity', 'P091308_rep1_2 MaxLFQ Intensity', 'P091308_rep2_3 MaxLFQ Intensity',
    'P091309_4 MaxLFQ Intensity', 'P091310_rep1_5 MaxLFQ Intensity', 'P091310_rep2_6 MaxLFQ Intensity',
    'P091311_7 MaxLFQ Intensity', 'P091312_8 MaxLFQ Intensity', 'P091313_9 MaxLFQ Intensity',
    'P091314_10 MaxLFQ Intensity', 'P091315_11 MaxLFQ Intensity', 'P091316_12 MaxLFQ Intensity',
    'P091317_13 MaxLFQ Intensity', 'P091318_14 MaxLFQ Intensity', 'P091319_rep1_15 MaxLFQ Intensity',
    'P091319_rep2_16 MaxLFQ Intensity', 'P091320_rep1_17 MaxLFQ Intensity', 'P091320_rep2_18 MaxLFQ Intensity',
    'P091321_19 MaxLFQ Intensity', 'P091322_20 MaxLFQ Intensity', 'P091323_21 MaxLFQ Intensity',
    'P091324_22 MaxLFQ Intensity', 'P091325_23 MaxLFQ Intensity', 'P091326_24 MaxLFQ Intensity',
    'P091327_25 MaxLFQ Intensity', 'P091328_26 MaxLFQ Intensity', 'P091329_27 MaxLFQ Intensity',
    'P091330_28 MaxLFQ Intensity', 'P094332_1 MaxLFQ Intensity', 'P094333_2 MaxLFQ Intensity',
    'P094334_3 MaxLFQ Intensity', 'P094335_4 MaxLFQ Intensity', 'P094336_5 MaxLFQ Intensity',
    'P094337_6 MaxLFQ Intensity', 'P094338_7 MaxLFQ Intensity', 'P094339_8 MaxLFQ Intensity',
    'P094340_9 MaxLFQ Intensity', 'P094341_10 MaxLFQ Intensity', 'P094342_11 MaxLFQ Intensity',
    'P094343_12 MaxLFQ Intensity', 'P094344_13 MaxLFQ Intensity', 'P094345_14 MaxLFQ Intensity',
    'P103997_1 MaxLFQ Intensity', 'P103998_2 MaxLFQ Intensity', 'P103999_3 MaxLFQ Intensity',
    'P104000_4 MaxLFQ Intensity', 'P104001_5 MaxLFQ Intensity', 'P104002_6 MaxLFQ Intensity',
    'P104003_7 MaxLFQ Intensity', 'P104004_8 MaxLFQ Intensity', 'P093501_1 MaxLFQ Intensity',
    'P093502_2 MaxLFQ Intensity', 'P093503_3 MaxLFQ Intensity', 'P093504_4 MaxLFQ Intensity',
    'P093505_5 MaxLFQ Intensity', 'P093506_6 MaxLFQ Intensity', 'P093507_7 MaxLFQ Intensity',
    'P093508_8 MaxLFQ Intensity', 'P093509_9 MaxLFQ Intensity', 'P093510_10 MaxLFQ Intensity',
    'P093511_11 MaxLFQ Intensity', 'P093512_12 MaxLFQ Intensity', 'P093513_13 MaxLFQ Intensity',
    'P093514_14 MaxLFQ Intensity', 'P093515_15 MaxLFQ Intensity', 'P093516_16 MaxLFQ Intensity']


# # Step 1: Mean of LFQ (tissue) values per HOG–Species group
# tissue_means = df.groupby(['HOG', 'Species'])[lfq_cols].mean().reset_index()

# # Step 2: Comma-separated protein IDs per group
# protein_lists = df.groupby(['HOG', 'Species'])['ProteinID'].apply(
#     lambda x: ','.join(pd.unique(x))
# ).reset_index()

# # Step 3: Merge both results
# collapsed_df = pd.merge(tissue_means, protein_lists, on=['HOG', 'Species'])

# # Optional: Rearranging columns
# collapsed_df = collapsed_df[['HOG', 'Species', 'ProteinID'] + lfq_cols]
# # Output to verify
# collapsed_df.to_csv("E:\Guido\sibel\Masters_Thesis\data\collapsed_firststep.tsv", sep='\t', index=False)
# print(collapsed_df)


import time, pandas as pd
from pathlib import Path

SOURCE = Path(r"E:\Guido\sibel\Masters_Thesis\data\hog_intensity_long.tsv")
OUT    = Path(r"E:\Guido\sibel\Masters_Thesis\data\collapsed_firststep.tsv")

t0 = time.perf_counter()
df = pd.read_csv(SOURCE, sep="\t", dtype=str)
lfq_cols = [c for c in df.columns if "MaxLFQ" in c]

# ---------- fast reducer that always casts to str --------------------
def join_str(series):
    return ",".join(series.dropna().astype(str))

list_aggs = {c: join_str for c in lfq_cols}

collapsed = (
    df.groupby(["HOG", "Species"], observed=True)
      .agg({"ProteinID": join_str, **list_aggs})
      .reset_index()
)

collapsed = collapsed[["HOG", "Species", "ProteinID"] + lfq_cols]
collapsed.to_csv(OUT, sep="\t", index=False)
print(f"✓ wrote {OUT}  ({time.perf_counter()-t0:.1f}s total)")


# Group by HOG and aggregate as required
collapsed = df.groupby('HOG').agg({
    'Species': lambda x: ','.join(sorted(pd.unique(x))),
    'ProteinID': lambda x: ','.join(pd.unique(','.join(x).split(','))),
    **{col: lambda x: ','.join(str(v) for v in x) for col in lfq_cols}
}).reset_index()
collapsed.to_csv("E:\Guido\sibel\Masters_Thesis\data\collapsed_final_HOG.tsv", sep='\t', index=False)