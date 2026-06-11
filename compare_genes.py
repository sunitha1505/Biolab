from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

genes = {
    "TP53": "data/tp53.fasta",
    "BRCA1": "data/brca1.fasta",
    "KRAS": "data/kras.fasta"
}

print("=" * 60)
print("GENE COMPARISON REPORT")
print("=" * 60)

for gene_name, filepath in genes.items():
    for record in SeqIO.parse(filepath, "fasta"):
        gc = round(gc_fraction(record.seq) * 100, 2)
        at = round(100 - gc, 2)
        print(f"\nGene: {gene_name}")
        print(f"Accession: {record.id}")
        print(f"Length: {len(record.seq)} bp")
        print(f"GC Content: {gc}%")
        print(f"AT Content: {at}%")

print("\n" + "=" * 60)

