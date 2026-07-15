import pubchempy as pcp
import os

# Neem phytochemicals + reference compound
compounds = {
    "Nimbolide":    12313376,
    "Azadirachtin": 5281303,
    "Gedunin":      12004512,
    "Nimbin":       108058,
    "Quercetin":    5280343
}

os.makedirs("data/compounds", exist_ok=True)

print("=" * 55)
print("NEEM PHYTOCHEMICAL FETCH — PubChem")
print("=" * 55)

results = []

for name, cid in compounds.items():
    print(f"\nFetching {name} (CID: {cid})...")
    
    compound = pcp.Compound.from_cid(cid)
    
    mw = compound.molecular_weight
    formula = compound.molecular_formula
    smiles = compound.smiles
    hbd = compound.h_bond_donor_count
    hba = compound.h_bond_acceptor_count
    logp = compound.xlogp
    rotatable = compound.rotatable_bond_count
    
    # Lipinski Rule of Five check
    violations = 0
    if mw > 500: violations += 1
    if logp > 5: violations += 1
    if hbd > 5: violations += 1
    if hba > 10: violations += 1
    lipinski = "PASS" if violations <= 1 else "FAIL"
    
    print(f"  Formula:    {formula}")
    print(f"  MW:         {mw} Da")
    print(f"  LogP:       {logp}")
    print(f"  HBD:        {hbd}")
    print(f"  HBA:        {hba}")
    print(f"  Rotatable:  {rotatable}")
    print(f"  Lipinski:   {lipinski}")
    
    # Save SMILES to file
    smiles_path = f"data/compounds/{name.lower()}.smi"
    with open(smiles_path, "w") as f:
        f.write(f"{smiles} {name}\n")
    print(f"  Saved to:   {smiles_path}")
    
    results.append({
        "name": name,
        "cid": cid,
        "formula": formula,
        "mw": mw,
        "logp": logp,
        "hbd": hbd,
        "hba": hba,
        "rotatable": rotatable,
        "lipinski": lipinski,
        "smiles": smiles
    })

print("\n" + "=" * 55)
print("SUMMARY TABLE")
print("=" * 55)
print(f"{'Compound':<15} {'MW':>8} {'LogP':>6} {'HBD':>5} {'HBA':>5} {'Lipinski':>10}")
print("-" * 55)
for r in results:
    print(f"{r['name']:<15} {r['mw']:>8} {r['logp']:>6} {r['hbd']:>5} {r['hba']:>5} {r['lipinski']:>10}")