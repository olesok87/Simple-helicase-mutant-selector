import pandas as pd
import glob
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

# 🔄 Load and sanitize each file
dfs = []
for i, file in enumerate(csv_files, start=1):
    df = pd.read_csv(file)

    # 🧼 Normalize column names
    df.columns = [col.strip().lower() for col in df.columns]

    # ✅ Check for required columns
    if not all(col in df.columns for col in shared_cols):
        raise ValueError(f"Missing one of the required columns in {file}: {shared_cols}")

    # 🔁 Convert resname to 1-letter code
    df['resname'] = df['resname'].str.upper().map(aa_map)
    if df['resname'].isnull().any():
        raise ValueError(f"Unrecognized amino acid code in {file}")

    # 🏷 Rename extra columns to keep them unique
    extra_cols = [col for col in df.columns if col not in shared_cols]
    df.rename(columns={col: f"{col}_f{i}" for col in extra_cols}, inplace=True)

    dfs.append(df)

# 🔗 Merge all DataFrames on shared columns
merged_df = reduce(lambda left, right: pd.merge(left, right, on=shared_cols, how='outer'), dfs)

# 💾 Save to CSV
merged_df.to_csv("merged_output.csv", index=False)
print("✅ Merged CSV saved as 'merged_output.csv'")