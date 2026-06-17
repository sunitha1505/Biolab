from Bio import SeqIO

# Load EGFR protein sequence
record = SeqIO.read("data/egfr_protein.fasta", "fasta")
protein = str(record.seq)

print("=" * 50)
print("EGFR MUTATION SITE ANALYSIS")
print("=" * 50)

# Position 858 — L858R mutation
pos_858 = protein[857] 
print(f"\nPosition 858: {pos_858}")
print(f"Normal: Leucine (L)")
print(f"Cancer: Arginine (R) — L858R mutation")
print(f"Effect: Locks kinase domain permanently ON")

# Position 790 — T790M resistance mutation
pos_790 = protein[789]
print(f"\nPosition 790: {pos_790}")
print(f"Normal: Threonine (T)")
print(f"Cancer: Methionine (M) — T790M mutation")
print(f"Effect: Narrows ATP pocket, blocks Erlotinib/Gefitinib")

print("\n" + "=" * 50)
print("Drug Response:")
print("L858R alone    → Erlotinib/Gefitinib works")
print("L858R + T790M  → Osimertinib needed (covalent)")
print("=" * 50)