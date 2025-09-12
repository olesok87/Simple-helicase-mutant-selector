#!/usr/bin/env python3
import csv
import sys
from Bio.PDB import MMCIFParser

def get_instability_residues(cif_path: str,
                             chain_id: str):
    """
    Scan a pre-filtered .exposed.cif for “problematic” residues
    in a single chain, based on the instability_map.

    cif_path : path to your .exposed.cif (only exposed residues present)
    chain_id : e.g. "A"
    """
    # map of “problematic” residues → (instability type, mutation suggestion)
    instability_map = {
        "ASN": ("Deamidation hotspot", "Mutate to D or A"),
        "GLN": ("Deamidation hotspot", "Mutate to E or A"),
        "MET": ("Oxidation-prone",      "Mutate to L or I"),
        "CYS": ("Free thiol",           "Mutate to S or A"),
        "PHE": ("Hydrophobic patch",    "Mutate to S or Y"),
        "TRP": ("Hydrophobic patch",    "Mutate to F or Y"),
        "TYR": ("Hydrophobic patch",    "Mutate to S or F"),
        "LEU": ("Hydrophobic patch",    "Mutate to T or S"),
        "ILE": ("Hydrophobic patch",    "Mutate to T or S"),
        "LYS": ("Protease site",        "Mutate to R")
    }

    # 1. Parse CIF and extract model[0]
    parser    = MMCIFParser(QUIET=True)
    structure = parser.get_structure("exposed", cif_path)
    model     = structure[0]

    inst_residues = []
    suggestions   = []

    # 2. Loop only over residues in the chosen chain
    if chain_id not in model:
        sys.exit(f"ERROR: Chain {chain_id} not found in {cif_path}")

    for res in model[chain_id]:
        # skip hetero-residues, keep only standard amino acids
        if res.id[0] != " ":
            continue

        resname = res.resname
        if resname not in instability_map:
            continue

        inst_type, suggestion = instability_map[resname]
        chain, resnum = chain_id, res.id[1]

        inst_residues.append(
            (chain, resnum, resname, inst_type)
        )
        suggestions.append(
            (chain, resnum, resname, suggestion)
        )

    return inst_residues, suggestions


### Helpers: write to CSV / TXT (unchanged)
def write_list_to_csv(filepath, data, headers):
    with open(filepath, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(headers)
        writer.writerows(data)
    print(f"→ Wrote {len(data)} records to {filepath}")

def write_suggestions_to_txt(filepath, suggestions):
    with open(filepath, "w") as fh:
        for chain, resnum, resname, suggestion in suggestions:
            fh.write(f"{chain} {resnum} {resname}: {suggestion}\n")
    print(f"→ Wrote {len(suggestions)} suggestions to {filepath}")


if __name__ == "__main__":
    cif_path = input("Path to .exposed.cif file: ").strip()
    chain_id = input("Chain to analyze (e.g. A): ").strip()

    inst_list, suggestion_list = get_instability_residues(cif_path, chain_id)

    write_list_to_csv(
        "/results/instability_residues.csv",
        inst_list,
        headers=["Chain","ResNum","ResName","Instability"]
    )
    write_suggestions_to_txt(
        "mutation_suggestions.txt",
        suggestion_list
    )