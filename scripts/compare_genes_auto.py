import sys
from Bio.SeqUtils import gc_fraction

sys.path.append("scripts")
from fetch_sequence import search_gene, fetch_record, extract_protein, analyze_protein 

genes = [gene.upper() for gene in sys.argv[1:]]

results = []

for gene in genes:
    print(f"Fetching {gene}...")

    ncbi_id = search_gene(gene)
    if ncbi_id is None:
        print(f" Skipping {gene} - not found")
        continue
    record = fetch_record(ncbi_id)

    gc = round(gc_fraction(record.seq) * 100, 2)
    at = round(100 - gc, 2)
    protein_seq = extract_protein(record)

    if protein_seq:
        mw, pi, instability, stability = analyze_protein(protein_seq)
        prot_length = len(protein_seq)
    else:
          mw, pi, instability, stability, prot_length = 0, 0, 0, "Unknown", 0
    
    results.append({
        "gene": gene,
        "length": len(record.seq),
        "gc": gc,
        "at": at,
        "protein_length": prot_length,
        "mw": round(mw/1000, 1),
        "pi": round(pi, 2),
        "instability": round(instability, 2),
        "stability": stability
    })

    # Step 4: print comparison table
print("\n" + "=" * 75)
print(f"{'Gene':<8} {'Length':>8} {'GC%':>7} {'Protein':>8} {'MW(kDa)':>9} {'pI':>6} {'Stability':<10}")
print("=" * 75)

for r in results:
    print(f"{r['gene']:<8} {r['length']:>8} {r['gc']:>7} {r['protein_length']:>8} {r['mw']:>9} {r['pi']:>6} {r['stability']:<10}")

# Step 5: save to file
with open("results/comparison_table.txt", "w") as f:
    f.write("GENE COMPARISON TABLE\n")
    f.write("=" * 75 + "\n")
    f.write(f"{'Gene':<8} {'Length':>8} {'GC%':>7} {'Protein':>8} {'MW(kDa)':>9} {'pI':>6} {'Stability':<10}\n")
    f.write("=" * 75 + "\n")
    for r in results:
        f.write(f"{r['gene']:<8} {r['length']:>8} {r['gc']:>7} {r['protein_length']:>8} {r['mw']:>9} {r['pi']:>6} {r['stability']:<10}\n")

print("\nSaved to results/comparison_table.txt")
          
