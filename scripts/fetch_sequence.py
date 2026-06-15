"""
fetch_sequence.py
Fetch any human gene sequence from NCBI by name — no manual downloads needed.

Usage:   python scripts/fetch_sequence.py <gene_name>
Example: python scripts/fetch_sequence.py TP53
         python scripts/fetch_sequence.py EGFR
"""

import sys
import os
import time
from Bio import Entrez, SeqIO

Entrez.email = "sunithathilak15@gmail.com"

def search_gene(gene_name, organism="Homo sapiens"):
    """Step 1: Search NCBI by gene name, return the top matching NCBI ID."""
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
    """Step 2: Download the full GenBank record for a given NCBI ID."""
    time.sleep(0.4)
    print(f"  Downloading record {ncbi_id} ...")

    handle = Entrez.efetch(db="nucleotide", id=ncbi_id, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()

    return record


def save_fasta(record, gene_name):
    """Step 3: Save sequence to data/<gene>_fetched.fasta"""
    os.makedirs("data", exist_ok=True)
    path = f"data/{gene_name.lower()}_fetched.fasta"

    with open(path, "w") as f:
        SeqIO.write(record, f, "fasta")

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

    print(f"\n  Accession:   {record.name}")
    print(f"  Description: {record.description}")
    print(f"  Length:      {len(record.seq)} bp")
    print(f"  Organism:    {record.annotations.get('organism', 'Unknown')}")

    path = save_fasta(record, gene)
    print(f"\n  Saved to: {path}")
    print("\nDone!")


if __name__ == "__main__":
    main()