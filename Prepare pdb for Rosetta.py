import os
from Bio import PDB

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

# 🔧 Set your paths here
input_pdb = r"C:\Users\aszyk\PycharmProjects\Simple_helicase_mutant_selector\pdb\6YYE.pdb"
output_all_chains_path = r"C:\Users\aszyk\PycharmProjects\Simple_helicase_mutant_selector\pdb\cleaned\cleaned_all.pdb"
output_single_chain_path = r"C:\Users\aszyk\PycharmProjects\Simple_helicase_mutant_selector\pdb\cleaned\single_chain.pdb"
chain_to_keep = input("Which chain would you like to analyse later in Rosetta?")  # Set to your desired chain


# 🚀 Run the cleaner
clean_and_renumber_pdb(input_pdb, output_all_chains_path, output_single_chain_path, chain_to_keep)