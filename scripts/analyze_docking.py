import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

results = {
    "Compound":  ["Quercetin", "Gedunin", "Nimbin", "Nimbolide"],
    "Affinity":  [-9.894, -5.27, -4.595, -4.40],
    "MW":        [302.23, 482.6, 540.6, 466.5],
    "LogP":      [1.5, 4.2, 2.3, 2.2],
    "Lipinski":  ["PASS", "PASS", "PASS", "PASS"],
    "Type":      ["Reference", "Neem", "Neem", "Neem"]
}

df = pd.DataFrame(results)
df = df.sort_values("Affinity")

print("=" * 60)
print("DOCKING RESULTS — Neem Phytochemicals vs COX-2")
print("Target: Human COX-2 (PDB: 5F19)")
print("Method: AutoDock Vina 1.2.7")
print("=" * 60)
print(df.to_string(index=False))

print(f"\nBest neem compound: {df[df['Type']=='Neem'].iloc[0]['Compound']}")
print(f"Best affinity: {df[df['Type']=='Neem'].iloc[0]['Affinity']} kcal/mol")
print(f"Reference (Quercetin): {df[df['Type']=='Reference'].iloc[0]['Affinity']} kcal/mol")

df.to_csv("results/docking_results.csv", index=False)

fig, axes = plt.subplots(1, 2, figsize=(12, 5))
fig.suptitle("Neem Phytochemicals vs COX-2 — Docking Results", fontweight='bold')

colors = ['#2196F3' if t == 'Reference' else '#4CAF50' for t in df['Type']]
axes[0].barh(df['Compound'], df['Affinity'], color=colors)
axes[0].set_xlabel('Binding Affinity (kcal/mol)')
axes[0].set_title('Binding Affinity (more negative = stronger)')
axes[0].axvline(x=0, color='black', linewidth=0.5)

axes[1].scatter(df['MW'], df['Affinity'], 
                c=colors, s=100, zorder=5)
for i, row in df.iterrows():
    axes[1].annotate(row['Compound'], 
                    (row['MW'], row['Affinity']),
                    textcoords="offset points", xytext=(5,5))
axes[1].set_xlabel('Molecular Weight (Da)')
axes[1].set_ylabel('Binding Affinity (kcal/mol)')
axes[1].set_title('MW vs Binding Affinity')

plt.tight_layout()
plt.savefig("results/docking_chart.png", dpi=150, bbox_inches='tight')
print("\nChart saved to results/docking_chart.png")
print("CSV saved to results/docking_results.csv")