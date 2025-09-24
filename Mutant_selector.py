import os
import csv
from Bio import PDB
from Bio.PDB import PDBParser, NeighborSearch, Selection
import numpy as np
import sys


# 0. Prep pdf file (can be skipped)
def clean_and_renumber_pdb(input_pdb, output_all_chains, output_single_chain, target_chain=None):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("protein", input_pdb)

    io = PDB.PDBIO()

    class SelectAllChains(PDB.Select):
        def accept_residue(self, residue):
            return residue.id[0] == " "  # Exclude heteroatoms and waters
        def accept_atom(self, atom):
            return atom.element != "H"  # Remove hydrogens

    class SelectSingleChain(PDB.Select):
        def accept_chain(self, chain):
            return chain.id == target_chain
        def accept_residue(self, residue):
            return residue.id[0] == " "
        def accept_atom(self, atom):
            return atom.element != "H"

    # Renumber all residues starting from 1 per chain
    for model in structure:
        for chain in model:
            for i, residue in enumerate(chain.get_residues(), start=1):
                residue.id = (" ", i, " ")

    # Save all chains
    io.set_structure(structure)
    io.save(output_all_chains, select=SelectAllChains())
    print(f"✅ Saved cleaned PDB with all chains to: {output_all_chains}")

    # Save only the specified chain
    io.set_structure(structure)
    io.save(output_single_chain, select=SelectSingleChain())
    print(f"✅ Saved cleaned PDB with chain {target_chain} to: {output_single_chain}")

# 1. Flexible Residue Detection via B-Factors
def get_flexible_residues(pdb_path, bfactor_threshold=30):
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure("structure", pdb_path)
    results = []
    for model in struct:
        for chain in model:
            for residue in chain:
                if not residue.has_id("CA"):
                    continue
                b_vals = [atom.get_bfactor() for atom in residue]
                avg_b = np.mean(b_vals)
                if avg_b > bfactor_threshold:
                    results.append((chain.id,
                                    residue.id[1],
                                    residue.resname,
                                    round(avg_b, 2)))
    return results


# 2.DNA-Contact Residue Detection
def get_dna_contact_residues(pdb_path, dna_chain_ids, contact_radius=5.0):
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure("complex", pdb_path)
    all_atoms = Selection.unfold_entities(struct, 'A')
    ns = NeighborSearch(all_atoms)

    dna_atoms = []
    for model in struct:
        for chain in model:
            if chain.id in dna_chain_ids:
                for residue in chain:
                    dna_atoms.extend(residue.get_atoms())

    contacts = set()
    for atom in dna_atoms:
        for residue in ns.search(atom.coord, contact_radius, level='R'):
            chain_id = residue.get_parent().id
            if chain_id in dna_chain_ids:
                continue
            resnum = residue.id[1]
            resname = residue.resname
            contacts.add((chain_id, resnum, resname))
    return list(contacts)


# 3. Nucleotide-Proximal Residue Detection
def get_nucleotide_proximal_residues(pdb_path, ligand_names, radius=6.0):
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure("structure", pdb_path)
    all_atoms = Selection.unfold_entities(struct, 'A')
    ns = NeighborSearch(all_atoms)

    lig_atoms = []
    for model in struct:
        for chain in model:
            for residue in chain:
                if residue.id[0] != 'H':
                    continue
                if residue.resname.strip().upper() in ligand_names:
                    lig_atoms.extend(residue.get_atoms())

    proximals = set()
    for atom in lig_atoms:
        for residue in ns.search(atom.coord, radius, level='R'):
            parent = residue.get_parent()
            if parent.id[0] != ' ':
                continue
            chain_id = parent.get_parent().id
            resnum = parent.id[1]
            resname = parent.resname
            proximals.add((chain_id, resnum, resname))
    return list(proximals)


# 4. Manual ATP-Site Residue Selection
def get_manual_site_residues(pdb_path, chain_id, residue_numbers):
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure("structure", pdb_path)
    manual = []
    for model in struct:
        for chain in model:
            if chain.id != chain_id:
                continue
            for num in residue_numbers:
                try:
                    res = chain[num]
                    manual.append((chain.id, res.id[1], res.resname))
                except KeyError:
                    print(f"Warning: residue {num} not found in chain {chain.id}")
    return manual


#from Bio import PDB
import sys

def get_instability_residues(pdb_path: str, chain_id: str):
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

    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_path)
    model = structure[0]

    inst_residues = []

    # Check if chain exists
    if chain_id not in [chain.id for chain in model]:
        sys.exit(f"ERROR: Chain {chain_id} not found in {pdb_path}")

    for res in model[chain_id]:
        if res.id[0] != " ":  # Skip heteroatoms and waters
            continue

        resname = res.get_resname()
        if resname not in instability_map:
            continue

        inst_type, suggestion = instability_map[resname]
        resnum = res.id[1]

        inst_residues.append((chain_id, resnum, resname, inst_type, suggestion))

    return inst_residues

### Helpers: write to CSV / TXT
def write_list_to_csv(filepath, data, headers):
    with open(filepath, 'w', newline='') as fh:
        writer = csv.writer(fh)
        writer.writerow(headers)
        writer.writerows(data)
    print(f"→ Wrote {len(data)} records to {filepath}")


# EXECUTION FLOW
if __name__ == "__main__":
    # pdb_file = input("🔍 PDB file path: ").strip()
    # fix
    input_pdb = r"C:\Users\aszyk\PycharmProjects\Simple_helicase_mutant_selector\pdb\6YYE.pdb"
    pdb_file = r"C:\Users\aszyk\PycharmProjects\Simple_helicase_mutant_selector\pdb\6YYE.pdb"
    #results_dir = input("📂 Results folder path: ").strip()
    results_dir=r"C:\Users\aszyk\PycharmProjects\Simple_helicase_mutant_selector\results"
    os.makedirs(results_dir, exist_ok=True)

    # 0. Prep pdf file
    output_all_chains_path = r"C:\Users\aszyk\PycharmProjects\Simple_helicase_mutant_selector\pdb\cleaned\cleaned_all.pdb"
    output_single_chain_path = r"C:\Users\aszyk\PycharmProjects\Simple_helicase_mutant_selector\pdb\cleaned\single_chain.pdb"
    chain_to_keep = input("Which chain would you like to analyse later in Rosetta?")  # Set to your desired chain
    clean_and_renumber_pdb(input_pdb, output_all_chains_path, output_single_chain_path, chain_to_keep)


    # 1. Flexible residues
    if input("✏️ Would you like to generate flexibility hotspots? (y/n): ").strip().lower() == 'y':
        bcut = float(input("📈 B-factor threshold (e.g., 30): ").strip())
        flex = get_flexible_residues(pdb_file, bcut)
        out1 = os.path.join(results_dir, "flexible_residues.csv")
        write_list_to_csv(out1, flex, ["chain", "resnum", "resname", "avg_bfactor"])

    # 2. DNA-contact residues
    if input("✏️ Would you like to generate DNA contact residues? (y/n): ").strip().lower() == 'y':
        dna_input = input("🧬 DNA chain IDs (comma-sep): ").strip()
        dna_chains = [c.strip() for c in dna_input.split(',') if c.strip()]
        dnac = get_dna_contact_residues(pdb_file, dna_chains)
        out2 = os.path.join(results_dir, "dna_contact_residues.csv")
        write_list_to_csv(out2, dnac, ["chain", "resnum", "resname"])

    # 3. Nucleotide-proximal residues
    lig_input = input("🧪 Ligands (comma-sep, e.g., ATP,ADP,ANP) or press enter: ").strip()
    ligs = [l.strip().upper() for l in lig_input.split(',') if l.strip()]
    prox = get_nucleotide_proximal_residues(pdb_file, ligs)
    out3 = os.path.join(results_dir, "nucleotide_proximal_residues.csv")
    write_list_to_csv(out3, prox, ["chain", "resnum", "resname"])

    # 4. Manual ATP-site residues
    if input("✏️ Manual site selection? (y/n): ").strip().lower() == 'y':
        cid = input("Chain ID for manual selection: ").strip()
        nums_input = input("Please use UniProt or similar for residue numbers (comma-sep): ").strip()
        nums = [int(x) for x in nums_input.split(',') if x.strip().isdigit()]
        manual = get_manual_site_residues(pdb_file, cid, nums)
        out4 = os.path.join(results_dir, "manual_site_residues.csv")
        write_list_to_csv(out4, manual, ["chain", "resnum", "resname"])

    # 5. Instability hotspots
    if input("✏️ Would you like to generate stability hotspots? (y/n): ").strip().lower() == 'y':
        #cif_path = input("Path to .exposed.pdf file: ").strip()
        cif_path = r'C:\Users\aszyk\PycharmProjects\Simple_helicase_mutant_selector\pdb\6YYE_exposed_surface.pdb'
        chain_id = input("Chain to analyze (e.g. A): ").strip()
        inst_list = get_instability_residues(cif_path, chain_id)
        out5_csv = os.path.join(results_dir, "instability_residues.csv")
        write_list_to_csv(out5_csv,inst_list,headers=["Chain","ResNum","ResName","Instability","Mutate to..."])


    # 6. Merging files
    # 7 . Generate Rosetta modelling file



