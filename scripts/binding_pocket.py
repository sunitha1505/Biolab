from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
import numpy as np

parser = PDBParser(QUIET=True)
structure = parser.get_structure("TP53", "data/tp53_alphafold.pdb")

print("=" * 55)
print("BINDING POCKET ANALYSIS — TP53")
print("=" * 55)

# Calculate solvent accessible surface area (SASA)
sr = ShrakeRupley()
sr.compute(structure, level="R") # R = residue level

# Hydrophobic amino acids - important for pocket lining
hydrophobic = ["ALA", "VAL", "ILE", "LEU", "MET", "PHE", "TRP", "PRO"]

# Collect buried hydrophobic residues in DNA-binding domain (92-291)
print("\n--- Buried Hydrophobic Residues in DNA-binding Domain ---")
print("(Low SASA = buried = potential pocket lining)")
print(f"\n{'Residue':<10} {'Number':<10} {'SASA':<10}")
print("-" * 30)

buried_residues = []
chain_a = structure[0]['A']

for residue in chain_a:
    res_num = residue.id[1]
    res_name = residue.resname

     # Only look at DNA-binding domain
    if res_num < 92 or res_num > 291:
        continue
    
    # Only hydrophobic residues
    if res_name not in hydrophobic:
        continue
    
    sasa = residue.sasa

      # Buried = SASA less than 20
    if sasa < 20:
        buried_residues.append((res_name, res_num, sasa))
        print(f"{res_name:<10} {res_num:<10} {sasa:.2f}")

print(f"\nTotal buried hydrophobic residues: {len(buried_residues)}")

# Find clusters — residues within 25 residues of each other
print("\n--- Potential Pocket Clusters ---")
print("(Groups of buried hydrophobic residues close in sequence)")

clusters = []
current_cluster = [buried_residues[0]]

for i in range(1, len(buried_residues)):
    if buried_residues[i][1] - buried_residues[i-1][1] <= 8:
        current_cluster.append(buried_residues[i])
    else:
        if len(current_cluster) >= 2:
            clusters.append(current_cluster)
        current_cluster = [buried_residues[i]]

if len(current_cluster) >= 3:
    clusters.append(current_cluster)

for i, cluster in enumerate(clusters):
    res_nums = [r[1] for r in cluster]
    print(f"\n  Cluster {i+1}: residues {min(res_nums)}-{max(res_nums)}")
    print(f"  Size: {len(cluster)} residues")
    print(f"  Residues: {', '.join([f'{r[0]}{r[1]}' for r in cluster])}")        




