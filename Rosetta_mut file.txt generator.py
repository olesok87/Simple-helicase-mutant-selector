import re
import pandas as pd

# 20 standard proteinogenic amino acids
AA_LIST = list("ACDEFGHIKLMNPQRSTVWY")

def parse_suggestions(raw):
    """
    Parse strings like "Mutate to E or A" or "E or A" into ['E','A'].
    Returns None if raw is empty/NaN.
    """
    if pd.isna(raw):
        return None
    s = re.sub(r'(?i)^mutate to', '', str(raw)).strip()
    if not s:
        return None
    # first try to pull single-letter tokens
    letters = re.findall(r'\b[A-Z]\b', s)
    if letters:
        return letters
    # fallback: split on commas/spaces/etc.
    parts = re.split(r'[,\s/;|]+', s)
    return [p.upper() for p in parts if p.upper() in AA_LIST] or None

def generate_rosetta_mutfile(
    merged_csv: str,
    chain: str,
    output_txt: str = r'C:\Users\aszyk\PycharmProjects\Simple_helicase_mutant_selector\results\mut_file.txt'
):
    """
    merged_csv:   path to your merged_output.csv
    chain:        the chain ID to process (e.g. "A" or "D")
    output_txt:   filename for the Rosetta .mut file
    """
    df = pd.read_csv(merged_csv)

    # detect columns
    inst_col = next(c for c in df.columns if 'instability' in c.lower())
    mut_col  = next(c for c in df.columns if 'mutate to' in c.lower())

    # filter to your chain
    df_chain = df[df['chain'] == chain].copy()
    if df_chain.empty:
        raise ValueError(f"No residues found for chain '{chain}'")

    designs = []
    for _, row in df_chain.sort_values('resnum').iterrows():
        wt = row['resname'].upper()
        rn = int(row['resnum'])

        if pd.notna(row[inst_col]):
            # only use the suggested mutants
            sugg = parse_suggestions(row[mut_col])
            # skip if nothing parseable
            targets = [aa for aa in (sugg or []) if aa != wt]
        else:
            # mutate to all other 19 amino acids
            targets = [aa for aa in AA_LIST if aa != wt]

        for aa in targets:
            designs.append((wt, rn, aa))

    # write out mut_file.txt
    with open(output_txt, 'w') as out:
        out.write(f"total {len(designs)}\n")
        for wt, rn, aa in designs:
            out.write("1\n")
            out.write(f"{wt} {rn} {aa}\n")

    print(f"✅ Wrote {len(designs)} mutations to '{output_txt}'")

if __name__ == "__main__":
    # === USER PARAMETERS ===
    merged_csv = r"C:\Users\aszyk\PycharmProjects\Simple_helicase_mutant_selector\results\merged_output.csv"  # your merged CSV
    chain_id = str(input("Which chain would you like to model?"))
    generate_rosetta_mutfile(merged_csv, chain_id)

