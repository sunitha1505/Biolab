from Bio import SeqIO, Entrez
from Bio.SeqUtils.ProtParam import ProteinAnalysis

Entrez.email = "ss15052007@gmail.com"

def get_protein(accession):
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()
    for feature in record.features:
        if feature.type == "CDS":
            cds_seq = feature.extract(record.seq)
            protein = cds_seq.translate(to_stop=True)
            return str(protein)
    return None

genes = {
    "TP53": "NM_000546.6",
    "KRAS": "NM_004985.5"
}

for gene_name, accession in genes.items():
    print(f"\nFetching {gene_name}...")
    protein_seq = get_protein(accession)

    analysis = ProteinAnalysis(protein_seq)

    print(f"\n{'='*60}")
    print(f"Gene: {gene_name}")
    print(f"Protein length: {len(protein_seq)} aa")
    print(f"Molecular weight: {round(analysis.molecular_weight(), 2)} Da")
    print(f"Isoelectric point (pI): {round(analysis.isoelectric_point(), 2)}")
    print(f"Instability index: {round(analysis.instability_index(), 2)}")
    stable = "Stable" if analysis.instability_index() < 40 else "Unstable"
    print(f"Stability: {stable}")