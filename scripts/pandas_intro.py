import pandas as pd

data = {
    "Gene":    ["TP53", "BRCA1", "KRAS", "EGFR"],
    "Length":  [2512, 7088, 5306, 9905],
    "GC":      [53.38, 41.77, 36.28, 47.78],
    "Protein": [393, 1863, 188, 1210],
    "MW_kDa":  [43.7, 207.7, 21.4, 134.3],
    "pI":      [6.33, 5.29, 8.24, 6.26],
    "Stability": ["Unstable", "Unstable", "Unstable", "Unstable"]
}

df = pd.DataFrame(data)

print("=== Full Table ===")
print(df)

print("\n=== Basic Statistics ===")
print(df.describe())

print("\n=== Sorted by GC Content ===")
print(df.sort_values("GC", ascending=False))

print("\n=== Genes with GC above 45% ===")
print(df[df["GC"] > 45])

print("\n=== Sorted by Protein Size ===")
print(df.sort_values("Protein", ascending=False)[["Gene", "Protein", "MW_kDa"]])

print("\n=== Average values ===")
print(f"Mean GC%: {df['GC'].mean():.2f}%")
print(f"Mean pI: {df['pI'].mean():.2f}")
print(f"Largest protein: {df.loc[df['Protein'].idxmax(), 'Gene']}")
print(f"Most GC-rich: {df.loc[df['GC'].idxmax(), 'Gene']}")

# Add mutation data
df['Key_mutations'] = ['R175H, R248W', 'BRCA1 185delAG', 'G12D, G12V', 'L858R, T790M']
print(df[['Gene', 'Key_mutations', 'MW_kDa']])

# Save to CSV
df.to_csv("results/gene_comparison.csv", index=False)
print("\nSaved to results/gene_comparison.csv")