from Bio.PDB import PDBParser
import numpy as np

# Load the structure
parser = PDBParser(QUIET=True)
structure = parser.get_structure("TP53", "data/tp53_structure.pdb")

# Basic info
print("=" * 50)
print("PDB STRUCTURE ANALYSIS — 3KMD")
print("TP53 Core Domain bound to DNA")
print("=" * 50)

# Loop through models, chains, residues
for model in structure:
    print(f"\nModel: {model.id}")
    for chain in model:
        residues = list(chain.get_residues())
        print(f"  Chain {chain.id}: {len(residues)} residues")

# Count atoms in each chain
print("\n--- Atom counts ---")
for model in structure:
    for chain in model:
        atoms = list(chain.get_atoms())
        print(f"  Chain {chain.id}: {len(atoms)} atoms")

# Look at first residue of Chain A
print("\n--- First residue of Chain A ---")
chain_a = structure[0]['A']
residues = list(chain_a.get_residues())
first_res = residues[0]
print(f"  Residue name: {first_res.resname}")
print(f"  Residue number: {first_res.id[1]}")
print(f"  Atoms in this residue:")
for atom in first_res:
    print(f"    {atom.name}: {atom.coord}")

# Find Arg248
print("\n--- Arg248 in Chain A (most mutated TP53 residue) ---")
chain_a = structure[0]['A']
arg248_nh2 = None
for residue in chain_a:
    if residue.id[1] == 248:
        print(f"  Residue: {residue.resname} {residue.id[1]}")
        for atom in residue:
            print(f"    {atom.name}: {atom.coord}")
            if atom.name == "NH2":
                arg248_nh2 = atom.coord

# Distance: Arg248 to DNA (both chains)
print("\n--- Distance: Arg248 to DNA (both chains) ---")
min_distance = 999
closest_atom = None
closest_residue = None
closest_chain = None

for chain_id in ['E', 'F']:
    dna_chain = structure[0][chain_id]
    for residue in dna_chain:
        if residue.resname == "HOH":
            continue
        for atom in residue:
            diff = arg248_nh2 - atom.coord
            distance = np.sqrt(np.sum(diff**2))
            if distance < min_distance:
                min_distance = distance
                closest_atom = atom.name
                closest_residue = residue.resname
                closest_chain = chain_id

print(f"  Arg248 NH2 position: {arg248_nh2}")
print(f"  Closest DNA atom: {closest_atom} in {closest_residue} (Chain {closest_chain})")
print(f"  Distance: {min_distance:.2f} Angstroms")

if min_distance < 3.5:
    print("  → Within hydrogen bonding distance (< 3.5 Å)")
elif min_distance < 5.0:
    print("  → Close contact but not direct hydrogen bond")
else:
    print("  → Too far for direct contact")