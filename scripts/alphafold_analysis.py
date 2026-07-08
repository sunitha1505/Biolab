from Bio.PDB import PDBParser
import numpy as np

parser = PDBParser(QUIET=True)

# Load AlphaFold structure
af_structure = parser.get_structure("TP53_AF", "data/tp53_alphafold.pdb")

print("=" * 55)
print("ALPHAFOLD vs CRYSTAL STRUCTURE — TP53")
print("=" * 55)

# Basic info
model = af_structure[0]
chain = model['A']
residues = [r for r in chain if r.resname != "HOH"]

print(f"\nAlphaFold Structure:")
print(f"  Total residues: {len(residues)}")
print(f"  Coverage: {residues[0].id[1]} to {residues[-1].id[1]}")

# pLDDT scores — stored in B-factor column in AlphaFold PDB files
print(f"\n--- pLDDT Confidence Scores ---")
plddt_scores = []
for residue in residues:
    for atom in residue:
        if atom.name == "CA":  # one score per residue via CA atom
            plddt_scores.append((residue.id[1], residue.resname, atom.bfactor))
            break

print(f"  Mean pLDDT: {np.mean([s[2] for s in plddt_scores]):.1f}")
print(f"  Min pLDDT:  {np.min([s[2] for s in plddt_scores]):.1f}")
print(f"  Max pLDDT:  {np.max([s[2] for s in plddt_scores]):.1f}")

# Count confidence categories
very_high = sum(1 for s in plddt_scores if s[2] > 90)
confident = sum(1 for s in plddt_scores if 70 < s[2] <= 90)
low = sum(1 for s in plddt_scores if 50 < s[2] <= 70)
very_low = sum(1 for s in plddt_scores if s[2] <= 50)

print(f"\n--- Confidence Distribution ---")
print(f"  Very high (>90):  {very_high} residues")
print(f"  Confident (70-90): {confident} residues")
print(f"  Low (50-70):       {low} residues")
print(f"  Very low (<50):    {very_low} residues")

# Check Arg248 confidence
print(f"\n--- Arg248 pLDDT (cancer hotspot) ---")
for res_num, resname, score in plddt_scores:
    if res_num == 248:
        print(f"  Residue: {resname} {res_num}")
        print(f"  pLDDT: {score}")
        if score > 90:
            print(f"  → Very high confidence — reliable for drug docking")
        elif score > 70:
            print(f"  → Confident — usable for analysis")
        else:
            print(f"  → Low confidence — treat with caution")

# Compare coverage with crystal structure
print(f"\n--- Comparison with Crystal Structure (3KMD) ---")
print(f"  Crystal structure: residues 92-291 (DNA-binding domain only)")
print(f"  AlphaFold:         residues {residues[0].id[1]}-{residues[-1].id[1]} (full protein)")
print(f"  AlphaFold covers {len(residues) - 200} more residues than crystal structure")