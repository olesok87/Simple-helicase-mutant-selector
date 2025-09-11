import csv
from Bio.PDB import PDBParser, NeighborSearch, Selection
import numpy as np

### 1️⃣ Flexible Residue Detection via B-Factors
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
                    results.append((chain.id, residue.id[1], residue.resname, round(avg_b, 2)))
    return results

# DNA-Contact Residue Detection (User-Defined Chains)
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
            resnum  = residue.id[1]
            resname = residue.resname
            contacts.add((chain_id, resnum, resname))
    return list(contacts)

# Nucleotide-Proximal Residue Detection (User-Defined Ligands)
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
            resnum   = parent.id[1]
            resname  = parent.resname
            proximals.add((chain_id, resnum, resname))
    return list(proximals)

# Manual ATP-Site Residue Selection
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

### Helper: write list of tuples to CSV
def write_list_to_csv(filename, data, headers):
    with open(filename, 'w', newline='') as fh:
        writer = csv.writer(fh)
        writer.writerow(headers)
        for row in data:
            writer.writerow(row)
    print(f"→ Wrote {len(data)} records to {filename}")

### ─────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    pdb_file = input("🔍 Enter path to your PDB file: ").strip()

    # 1. Flexible residues
    bcut = float(input("📈 B-factor threshold (e.g., 30): ").strip())
    flex = get_flexible_residues(pdb_file, bcut)
    print(f"\n🧬 Flexible residues found: {len(flex)}")
    write_list_to_csv("flexible_residues.csv", flex,
                      ["chain", "resnum", "resname", "avg_bfactor"])

    # 2. DNA-contact residues
    dna_input = input("\n🧬 DNA chain IDs (comma-sep, e.g., D,E): ").strip()
    dna_chains = [c.strip() for c in dna_input.split(',') if c.strip()]
    dna_contacts = get_dna_contact_residues(pdb_file, dna_chains)
    print(f"🧲 DNA-contact residues found: {len(dna_contacts)}")
    write_list_to_csv("dna_contact_residues.csv", dna_contacts,
                      ["chain", "resnum", "resname"])

    # 3. Nucleotide-proximal residues
    lig_input = input("\n🧪 Nucleotide ligands (comma-sep, e.g., ATP,ADP,ANP): Press enter if none ").strip()
    ligs = [l.strip().upper() for l in lig_input.split(',') if l.strip()]
    prox = get_nucleotide_proximal_residues(pdb_file, ligs)
    print(f"⚡ Nucleotide-proximal residues found: {len(prox)}")
    write_list_to_csv("nucleotide_proximal_residues.csv", prox,
                      ["chain", "resnum", "resname"])

    # 4. Manual site selection
    manual_q = input("\n✏️ Manually select residues around ATP site? (y/n): ").strip().lower()
    if manual_q == 'y':
        chain_id = input("🔢 Enter chain ID for manual selection (e.g., A): ").strip()
        nums_input = input("📋 Enter residue numbers based on UNIPROT, make sure numbering is the same (comma-sep, e.g., 1,2,5,8): ").strip()
        nums = []
        for x in nums_input.split(','):
            try:
                nums.append(int(x.strip()))
            except ValueError:
                print(f"Warning: '{x}' is not a valid integer—skipped.")
        manual = get_manual_site_residues(pdb_file, chain_id, nums)
        print(f"✂️ Manually selected residues: {len(manual)}")
        write_list_to_csv("manual_site_residues.csv", manual,
                          ["chain", "resnum", "resname"])