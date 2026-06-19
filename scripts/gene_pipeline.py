import sys
import pandas as pd
from datetime import datetime

sys.path.append("scripts")
from fetch_sequence import search_gene, fetch_record, extract_protein, analyze_protein
from Bio.SeqUtils import gc_fraction

def main():
    if len(sys.argv) < 2:
         print("Usage: python3 scripts/gene_pipeline.py <gene1> <gene2> ...")
         print("Example: python3 scripts/gene_pipeline.py TP53 BRCA1 KRAS EGFR")
         sys.exit(1)

    genes = [g.upper() for g in sys.argv[1:]]
    organism = "Homo sapiens"
    print(f"\nGene Pipeline — {organism}")
    print(f"Genes: {', '.join(genes)}")
    print("=" * 50)

# Step 2 & 3: fetch each gene and store results
    results = []

    for gene in genes:
        print(f"\nFetching {gene}...")
        
        ncbi_id = search_gene(gene, organism)
        if ncbi_id is None:
            print(f"  Skipping {gene} - not found")
            continue

        record = fetch_record(ncbi_id)
        gc = round(gc_fraction(record.seq) * 100, 2)
        at = round(100 - gc, 2)
        protein_seq = extract_protein(record)

        if protein_seq:
            mw, pi, instability, stability = analyze_protein(protein_seq)
            prot_len = len(protein_seq)
        else:
            mw, pi, instability, stability, prot_len = 0, 0, 0, "Unknown", 0

        results.append({
            "Gene":        gene,
            "Accession":   record.name,
            "mRNA_Length": len(record.seq),
            "GC_percent":  gc,
            "AT_percent":  at,
            "Protein_aa":  prot_len,
            "MW_kDa":      round(mw / 1000, 1),
            "pI":          round(pi, 2),
            "Instability": round(instability, 2),
            "Stability":   stability
        })
        print(f"  Done — {len(record.seq)} bp, {prot_len} aa, GC {gc}%")

    print("\nAll genes fetched.")

# Step 4: load into pandas
    df = pd.DataFrame(results)

    # Step 5: print table and statistics
    print("\n" + "=" * 70)
    print("RESULTS TABLE")
    print("=" * 70)
    print(df.to_string(index=False))

    print("\n" + "=" * 70)
    print("STATISTICS")
    print("=" * 70)
    print(f"Genes analyzed:     {len(df)}")
    print(f"Mean GC%:           {df['GC_percent'].mean():.2f}%")
    print(f"Mean protein size:  {df['Protein_aa'].mean():.0f} aa")
    print(f"Most GC-rich:       {df.loc[df['GC_percent'].idxmax(), 'Gene']}")
    print(f"Largest protein:    {df.loc[df['Protein_aa'].idxmax(), 'Gene']}")
    print(f"Highest pI:         {df.loc[df['pI'].idxmax(), 'Gene']} ({df['pI'].max()})")

    # Step 6: export to CSV with timestamp
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    csv_path = f"results/gene_pipeline_{timestamp}.csv"
    df.to_csv(csv_path, index=False)
    print(f"\nCSV saved to: {csv_path}")

    # Step 7: save summary txt
    txt_path = f"results/gene_pipeline_{timestamp}_summary.txt"
    with open(txt_path, "w") as f:
        f.write("GENE PIPELINE REPORT\n")
        f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M')}\n")
        f.write(f"Organism: {organism}\n")
        f.write(f"Genes: {', '.join(genes)}\n")
        f.write("=" * 70 + "\n")
        f.write(df.to_string(index=False))
        f.write("\n\n")
        f.write(f"Most GC-rich: {df.loc[df['GC_percent'].idxmax(), 'Gene']}\n")
        f.write(f"Largest protein: {df.loc[df['Protein_aa'].idxmax(), 'Gene']}\n")
    print(f"Summary saved to: {txt_path}")

main()
