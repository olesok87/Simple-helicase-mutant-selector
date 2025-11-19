Run
1. Prepare pdb
2. Select surface residues in pymol (please see below pymol commands) and dsave as pdb
3. Mutant_selector.py
4. Merging.py
5. Rosetta_mut file.txt generator.py

Pymol surface exposure commands:
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




