import freesasa
from Bio.PDB import PDBParser
import os
from Bio.PDB.PDBExceptions import PDBConstructionException


# --- Configuration ---
# Set the path to your PDB file. Make sure this path is correct.
pdb_file_path = r"C:\Users\aszyk\PycharmProjects\Simple_helicase_mutant_selector\pdb\6Y6C.pdb"
struct = freesasa.structureFromBioPDB(pdb_file_path, classifier={})

# Define the SASA (Solvent Accessible Surface Area) threshold in Å²
# to consider an amino acid as a "surface" residue.
# Residues with SASA > sasa_threshold will be identified as surface amino acids.
# Note: The provided sources do not define what constitutes a "surface amino acid"
# or suggest a specific SASA threshold. This value is user-defined based on
# common practices in structural bioinformatics or your specific research criteria.
sasa_threshold = 20.0

# --- Main Logic for Surface Amino Acid Identification ---

if not os.path.exists(pdb_file_path):
    print(f"Error: The PDB file '{pdb_file_path}' was not found.")
else:
    print(f"Loading and processing {pdb_file_path}...")

    try:
        # 1. Parse the PDB file using Bio.PDB.PDBParser.
        parser = PDBParser()

        # --- CRITICAL CORRECTION APPLIED HERE (THIRD ATTEMPT) ---
        # Extract a simple string ID for the structure (e.g., "6Y6C" from "6Y6C.pdb").
        # os.path.splitext() returns a tuple (root, ext), we take only the first element .
        structure_id = os.path.splitext(os.path.basename(pdb_file_path))  # This gets ('6Y6C', '.pdb')
        # We need the first element of that tuple, which is '6Y6C'
        structure_id = structure_id  # Now structure_id is '6Y6C'
        # ---------------------------------------------------------

        print("Step 1: Parsing PDB file with Bio.PDB using ID:", structure_id, "...")
        # As per the sources, parser.get_structure expects a string ID [1].
        bio_structure = parser.get_structure(structure_id, pdb_file_path)
        print("Step 1 Complete: Bio.PDB structure successfully parsed (warnings may appear above).")

        # 2. Convert the Bio.PDB structure to a freesasa.Structure object.
        # This provides more control over the analysis for freesasa calculations, as stated in the sources [2].
        print("Step 2: Converting Bio.PDB structure to freesasa.Structure...")
        freesasa_structure = freesasa.structureFromBioPDB(bio_structure)
        print("Step 2 Complete: freesasa.Structure successfully created.")

        # 3. Calculate the overall SASA for the entire protein.
        # This is a basic calculation using defaults, as shown in the sources [3].
        print("Step 3: Calculating overall SASA using freesasa.calc()...")
        result = freesasa.calc(freesasa_structure)
        print(f"Step 3 Complete: Overall SASA calculated. Total Area: {result.totalArea():.2f} A2")

        # 4. Prepare a list of selections for each individual standard amino acid residue.
        selections_to_calculate = []
        residue_map_for_print = {}

        print("Step 4: Preparing selections for individual amino acids...")
        for model in bio_structure:
            for chain in model:
                for residue in chain:
                    # Filter out heteroatoms (e.g., HOH for water, ions)
                    # Standard amino acid residues have a ' ' (space) as their hetero flag.
                    het_flag, res_seq, i_code = residue.get_id()
                    if het_flag != ' ':  # Skip non-standard residues (HETATM)
                        continue

                    res_name = residue.get_resname()
                    chain_id = chain.get_id()

                    # Create a unique label for this residue and a Pymol-like selection string.
                    # freesasa.selectArea() can use a subset of Pymol selection syntax [4].
                    label = f"{chain_id}_{res_name}_{res_seq}_{i_code}".strip('_')  # Clean up if i_code is empty
                    selection_string = f"resn {res_name} and resi {res_seq} and chain {chain_id}"

                    selections_to_calculate.append((label, selection_string))
                    residue_map_for_print[label] = (chain_id, res_name, res_seq, i_code)

        print(f"Step 4 Complete: Prepared {len(selections_to_calculate)} residue selections.")

        # 5. Calculate the SASA for each prepared residue selection.
        sasa_per_residue = {}
        if selections_to_calculate:
            print("Step 5: Calculating SASA for each residue using freesasa.selectArea()...")
            # freesasa.selectArea() integrates SASA over specified atom selections [4].
            sasa_per_residue = freesasa.selectArea(selections_to_calculate, freesasa_structure, result)
            print("Step 5 Complete: Individual residue SASA calculated.")
        else:
            print("No standard amino acid residues found in the PDB file to calculate SASA. Skipping Step 5.")

        # 6. Identify and print the surface amino acids based on the SASA threshold.
        print(f"\nIdentified Surface Amino Acids (SASA > {sasa_threshold:.2f} Å²):")
        print(f"{'Chain':<6} {'Residue':<10} {'ResID':<8} {'iCode':<7} {'SASA (A^2)':<15}")
        print("-" * 50)

        found_surface_aa = False
        for label, sasa_value in sasa_per_residue.items():
            if sasa_value > sasa_threshold:
                chain_id, res_name, res_seq, i_code = residue_map_for_print[label]

                print(f"{chain_id:<6} {res_name:<10} {str(res_seq):<8} {str(i_code):<7} {sasa_value:<15.2f}")
                found_surface_aa = True

        if not found_surface_aa:
            print("No surface amino acids found based on the specified threshold.")

    except PDBConstructionException as e:
        print(f"An error occurred during Bio.PDB parsing: {e}")
    except FileNotFoundError:
        print(f"Error: The file '{pdb_file_path}' was not found.")
    except Exception as e:
        print(f"\n--- CRITICAL ERROR ---")
        print(f"An unexpected error occurred during processing after Bio.PDB parsing. ")
        print(f"This error likely originated from the freesasa library itself.")
        print(f"Error details: {e}")
        print(f"Please review the PDB file '{pdb_file_path}' for any unusual features or severe discontinuities.")
        print(f"The warnings from Bio.PDB above might be related to the cause.")