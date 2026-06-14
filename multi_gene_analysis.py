from Bio import SeqIO, Entrez
from Bio.SeqUtils import gc_fraction
from Bio.SeqUtils.ProtParam import ProteinAnalysis

Entrez.email = "ss15052007@gmail.com"

def get_protein(accession):
    handle = Entrez.efetch(db="nucleotide", id=accession, 
                          rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()
    for feature in record.features:
        if feature.type == "CDS":
            cds_seq = feature.extract(record.seq)
            protein = cds_seq.translate(to_stop=True)
            return str(protein)
    return None

def analyze_nucleotide(fasta_path):
    for record in SeqIO.parse(fasta_path, "fasta"):
        gc = round(gc_fraction(record.seq) * 100, 2)
        at = round(100 - gc, 2)
        return {
            "length": len(record.seq),
            "gc": gc,
            "at": at
        }

def analyze_protein(protein_seq):
    analysis = ProteinAnalysis(protein_seq)
    stable = "Stable" if analysis.instability_index() < 40 else "Unstable"
    return {
        "length": len(protein_seq),
        "mw": round(analysis.molecular_weight(), 2),
        "pi": round(analysis.isoelectric_point(), 2),
        "instability": round(analysis.instability_index(), 2),
        "stability": stable
    }

genes = {
    "TP53":  ("data/tp53.fasta",  "NM_000546.6"),
    "BRCA1": ("data/brca1.fasta", "NM_007294.4"),
    "KRAS":  ("data/kras.fasta",  "NM_004985.5")
}

output = open("results/multi_gene_report.txt", "w")
output.write("=" * 70 + "\n")
output.write("MULTI-GENE ANALYSIS REPORT\n")
output.write("Genes: TP53, BRCA1, KRAS\n")
output.write("=" * 70 + "\n")

for gene_name, (fasta_path, accession) in genes.items():
    print(f"Analyzing {gene_name}...")

    nuc = analyze_nucleotide(fasta_path)
    protein_seq = get_protein(accession)
    prot = analyze_protein(protein_seq)

    output.write(f"\nGene: {gene_name}\n")
    output.write("-" * 40 + "\n")
    output.write(f"mRNA Length:        {nuc['length']} bp\n")
    output.write(f"GC Content:         {nuc['gc']}%\n")
    output.write(f"AT Content:         {nuc['at']}%\n")
    output.write(f"Protein Length:     {prot['length']} aa\n")
    output.write(f"Molecular Weight:   {prot['mw']} Da\n")
    output.write(f"Isoelectric Point:  {prot['pi']}\n")
    output.write(f"Instability Index:  {prot['instability']}\n")
    output.write(f"Stability:          {prot['stability']}\n")

output.write("\n" + "=" * 70 + "\n")
output.write("END OF REPORT\n")
output.close()

print("\nDone. Report saved to results/multi_gene_report.txt")