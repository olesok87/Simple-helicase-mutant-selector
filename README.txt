This script were designed to generate a list of mutants for a helicase (but can be used for any other protein):

Prepare pdb for Rosetta.py:
- prepare pdb for rosetta (cleaning, renumbering,atom renumbering)
- select the chain to be analysed by Rosetta to be used as input alongside with mut file (Rosetta ddG monomer protocol will process only one chain))



Mutant selector.py:
The primary residues to be mutated are:
- DNA contact residues (automatically extracted from .pdb after DNA chain is selected by user)
- Flexible residues (user selects b factor cutoff,automatically extracted from the pdb file)
- Instability residues (manually provided checking common instability residues) again surface residues calculated by pymol. Extra pdb will have to be provided by the user from pymol)
- Nucleotide proximal residues (automatic or manual)
- active site residues (manually provided by user)

Merging.py
- if needed this script is going to merge all mutants and remove redundancies to generate a final list

Rosetta_mut file.txt generator.py
- final list of mutants is converted to a rosetta ddg monomer mut file format

Run
1. Prepare pdb
2. Select surface residues in pymol (please see below pymol commands) and dsave as pdb
3. Mutant_selector.py
4. Merging.py
5. Rosetta_mut file.txt generator.py

* Pymol surface exposure commands:
Downlaod separate mmcif file (pdb didn't work)
# Compute SASA and store under B-factor
set dot_solvent, 1 ~(result " Setting: dot_solvent set to on.") solvent accsessibility is turned on
get_area all, load_b=1 # calculates solvent area and stores it in the b-factor column
#select for solvent exposed residues
select exposed_atoms, b > 30 (you can chnage the 30A^2 value here)
select exposed_residues, byres exposed_atoms (creates a residue selection from the atom selection)
color yellow, exposed_residues (visial help)
show sticks, exposed_residues
save exposed_surface.cif, exposed_residues # or export to object and save


ROSETTA ddG MONOMER PROTOCOL:
Rosetta flags for ddG monomer:

The structure was fixed with PDB Fixer to add missing atoms (web server) and Rosetta pre-minimised (no side chain movements due to size to start with) and later analysed by high resolution Rosetta protocol (sc are back on)
https://docs.rosettacommons.org/docs/latest/application_documentation/analysis/ddg-monomer


Rosetta pre-minimise flags
 /home/aszyk/rosetta/rosetta.binary.ubuntu.release-371/main/source/bin/minimize_with_cst.static.linuxgccrelease -in:file:s /home/aszyk/rosetta/helicase_inputs/2P6R_single_chain_PDB_Fixer_Repair.pdb -in:file:fullatom -ignore_unrecognized_res -fa_max_dis 9.0 -database /home/aszyk/rosetta/rosetta.binary.ubuntu.release-371/main/database/ -ddg::harmonic_ca_tether 0.5 -ddg::constraint_weight 1.0 -ddg::out_pdb_prefix min_cst_no_sc_ -ddg::sc_min_only true


I had to decrease some sampling here to increase the speed since my laptop is very slow. I would run it very differently on a cluster with more iterations.

/home/aszyk/rosetta/rosetta.binary.ubuntu.release-371/main/source/bin/ddg_monomer.static.linuxgccrelease \
    -s /home/aszyk/rosetta/helicase_inputs/output/min_cst_no_sc_.2P6R_single_chain_PDB_Fixer_Repair_0001.pdb \
    -ddg:mut_file /home/aszyk/rosetta/helicase_inputs/mut_file2.txt \
    -ddg:weight_file soft_rep_design \
    -fa_max_dis 6.0 # optional  \
    -ddg::iterations 20 \
    -ddg::dump_pdbs false \
    -ignore_unrecognized_res \
    -ddg::local_opt_only true \
    -ddg::min_cst true \
    -constraints: /home/aszyk/rosetta/helicase_inputs/output/helicase_constraints_new.cst \
    -in::file::fullatom \
    -ddg::mean false \
    -ddg::min true \
    -ddg::sc_min_only false \
    -ddg::ramp_repulsive true




