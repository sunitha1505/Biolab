"""
fetch_sequence.py
Fetch any human gene sequence from NCBI by name.

Usage:   python scripts/fetch_sequence.py <gene_name>
Example: python scripts/fetch_sequence.py TP53
"""

import sys
import os
import time
from Bio import Entrez, SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils import gc_fraction

Entrez.email = "your.email@example.com"


def search_gene(gene_name, organism="Homo sapiens"):
    """Search NCBI by gene name, return top NCBI ID."""
    query = f"{gene_name}[Gene Name] AND {organism}[Organism] AND mRNA[Filter] AND srcdb_refseq[Properties] AND refseq_select[Filter]"
    print(f"  Searching NCBI: {query}")

    handle = Entrez.esearch(db="nucleotide", term=query, retmax=1)
    result = Entrez.read(handle)
    handle.close()

    ids = result["IdList"]
    if not ids:
        print(f"  No results found for '{gene_name}'.")
        return None

    print(f"  Found NCBI ID: {ids[0]}")
    return ids[0]


def fetch_record(ncbi_id):
    """Download full GenBank record for a given NCBI ID."""
    time.sleep(0.4)
    print(f"  Downloading record {ncbi_id} ...")

    handle = Entrez.efetch(db="nucleotide", id=ncbi_id, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()

    return record


def save_fasta(record, gene_name):
    """Save nucleotide sequence to data/<gene>_fetched.fasta"""
    os.makedirs("data", exist_ok=True)
    path = f"data/{gene_name.lower()}_fetched.fasta"
    with open(path, "w") as f:
        SeqIO.write(record, f, "fasta")
    return path


def extract_protein(record):
    """Extract protein sequence from CDS features."""
    for feature in record.features:
        if feature.type == "CDS":
            if "translation" in feature.qualifiers:
                return feature.qualifiers["translation"][0]
    return None


def analyze_protein(protein_seq):
    """Compute protein properties."""
    pa = ProteinAnalysis(protein_seq)
    mw = pa.molecular_weight()
    pi = pa.isoelectric_point()
    instability = pa.instability_index()
    stability = "Stable" if instability < 40 else "Unstable"
    return mw, pi, instability, stability


def save_protein_fasta(protein_seq, gene_name, accession):
    """Save protein sequence to data/<gene>_protein.fasta"""
    path = f"data/{gene_name.lower()}_protein.fasta"
    with open(path, "w") as f:
        f.write(f">{accession} {gene_name} protein\n")
        f.write(f"{protein_seq}\n")
    return path


def main():
    if len(sys.argv) < 2:
        print("Usage: python scripts/fetch_sequence.py <gene_name>")
        print("Example: python scripts/fetch_sequence.py TP53")
        sys.exit(1)

    gene = sys.argv[1].upper()
    print(f"\nFetching {gene} from NCBI...")
    print("=" * 50)

    ncbi_id = search_gene(gene)
    if ncbi_id is None:
        sys.exit(1)

    record = fetch_record(ncbi_id)

    # --- Nucleotide ---
    gc = gc_fraction(record.seq) * 100
    print(f"\n--- Nucleotide ---")
    print(f"  Accession:   {record.name}")
    print(f"  Description: {record.description}")
    print(f"  Length:      {len(record.seq)} bp")
    print(f"  GC Content:  {gc:.2f}%")
    print(f"  Organism:    {record.annotations.get('organism', 'Unknown')}")
    nuc_path = save_fasta(record, gene)
    print(f"  Saved to:    {nuc_path}")

    # --- Protein ---
    protein_seq = extract_protein(record)
    if protein_seq:
        mw, pi, instability, stability = analyze_protein(protein_seq)
        print(f"\n--- Protein ---")
        print(f"  Length:      {len(protein_seq)} aa")
        print(f"  Mol. Weight: {mw/1000:.1f} kDa")
        print(f"  pI:          {pi:.2f}")
        print(f"  Instability: {instability:.2f} ({stability})")
        prot_path = save_protein_fasta(protein_seq, gene, record.name)
        print(f"  Saved to:    {prot_path}")
    else:
        print(f"\n  No protein translation found in record.")

    print("\nDone!")


if __name__ == "__main__":
    main()