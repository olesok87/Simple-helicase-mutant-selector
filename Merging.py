import os
import pandas as pd
from functools import reduce

# 🔧 Step 1: Define your input files
csv_files = [
    r"C:\Users\aszyk\PycharmProjects\Simple_helicase_mutant_selector\results\dna_contact_residues.csv",
    r"C:\Users\aszyk\PycharmProjects\Simple_helicase_mutant_selector\results\flexible_residues.csv",
    r"C:\Users\aszyk\PycharmProjects\Simple_helicase_mutant_selector\results\instability_residues.csv",
    r"C:\Users\aszyk\PycharmProjects\Simple_helicase_mutant_selector\results\manual_site_residues.csv",
    r"C:\Users\aszyk\PycharmProjects\Simple_helicase_mutant_selector\results\nucleotide_proximal_residues.csv"
]

# 🧬 3-letter to 1-letter amino acid mapping
aa_map = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}

# 🧹 Shared columns expected in all files
shared_cols = ['chain', 'resnum', 'resname']

# 🔄 Load, sanitize, and collect valid DataFrames
loaded_dfs = []
missing_files = []
invalid_files = []

for idx, path in enumerate(csv_files, start=1):
    if not os.path.exists(path):
        missing_files.append(path)
        continue

    try:
        df = pd.read_csv(path)
    except Exception as e:
        invalid_files.append((path, str(e)))
        continue

    # 🧼 Normalize column names
    df.columns = [col.strip().lower() for col in df.columns]

    # ✅ Ensure required columns are present
    if not all(col in df.columns for col in shared_cols):
        print(f"⚠️  Skipping {path}: missing one of {shared_cols}")
        continue

    # 🔁 Convert 3-letter resname to 1-letter code
    df['resname'] = df['resname'].str.upper().map(aa_map)
    if df['resname'].isnull().any():
        raise ValueError(f"Unrecognized amino acid code found in {path}")

    # 🏷 Prefix other columns so they stay unique post-merge
    extra_cols = [c for c in df.columns if c not in shared_cols]
    rename_map = {c: f"{c}_f{idx}" for c in extra_cols}
    df.rename(columns=rename_map, inplace=True)

    # annotate source and collect
    df['source_file'] = os.path.basename(path)
    loaded_dfs.append(df)

# 🚫 Report any files that didn’t load
if missing_files:
    print("\n🚫 The following files were not found and were skipped:")
    for f in missing_files:
        print(f"  - {f}")

if invalid_files:
    print("\n🚫 The following files failed to load and were skipped:")
    for f, err in invalid_files:
        print(f"  - {f}  (error: {err})")

# 🛑 If nothing loaded, stop here
if not loaded_dfs:
    raise RuntimeError("No valid CSV files available to merge.")

# 🔗 Merge all loaded DataFrames on shared columns
merged_df = reduce(
    lambda left, right: pd.merge(left, right, on=shared_cols, how='outer'),
    loaded_dfs
)

# 💾 Save the merged result
out_path = r"C:\Users\aszyk\PycharmProjects\Simple_helicase_mutant_selector\results\merged_output.csv"
merged_df.to_csv(out_path, index=False)
print(f"\n✅ Merged CSV saved as '{out_path}'.")