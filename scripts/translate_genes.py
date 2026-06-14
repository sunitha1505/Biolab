from Bio import Entrez, SeqIO
from Bio.SeqUtils import gc_fraction

Entrez.email = "sunithathilak15@email.com"

def fetch_gene_info(gene_name, accession):
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()

    gc = gc_fraction(record.seq) * 100

    # ✅ Find the CDS feature — don't translate the raw mRNA
    protein_seq = None
    for feature in record.features:
        if feature.type == "CDS":
            cds_seq = feature.extract(record.seq)   # extracts only the coding region
            protein_seq = cds_seq.translate(to_stop=True)
            break

    # Alternatively, use the pre-translated qualifier (faster, no warning):
    # protein_seq = feature.qualifiers["translation"][0]

    print("=" * 60)
    print(f"Gene: {gene_name}")
    print(f"Accession: {accession}")
    print(f"mRNA Length: {len(record.seq)} bp")
    print(f"GC Content: {gc:.2f}%")
    print(f"Protein Length: {len(protein_seq)} amino acids")
    print(f"First 10 amino acids: {protein_seq[:10]}")
    print(f"Starts with Met (M): {str(protein_seq).startswith('M')}")
    print("=" * 60)

fetch_gene_info("TP53", "NM_000546.6")
fetch_gene_info("KRAS", "NM_004985.5")