import os
import subprocess

compounds_dir = "data/compounds"
ligands_dir = "data/ligands"
os.makedirs(ligands_dir, exist_ok=True)

smi_files = [f for f in os.listdir(compounds_dir) if f.endswith(".smi")]

print("=" * 55)
print("PREPARING LIGANDS FOR DOCKING")
print("=" * 55)

for smi_file in smi_files:
    name = smi_file.replace(".smi", "")
    smi_path = os.path.join(compounds_dir, smi_file)
    sdf_path = os.path.join(ligands_dir, f"{name}.sdf")
    pdbqt_path = os.path.join(ligands_dir, f"{name}.pdbqt")

    print(f"\nProcessing {name}...")

    cmd1 = ["obabel", smi_path, "-O", sdf_path, "--gen3d", "fast", "-h"]
    try:
        result1 = subprocess.run(cmd1, capture_output=True, text=True, timeout=60)
    except subprocess.TimeoutExpired:
        print(f"  TIMEOUT — skipping {name}")
        continue

    if os.path.exists(sdf_path):
        print(f"  SDF generated: {sdf_path}")
    else:
        print(f"  ERROR generating SDF: {result1.stderr}")
        continue

    cmd2 = ["obabel", sdf_path, "-O", pdbqt_path, "-xh"]
    try:
        result2 = subprocess.run(cmd2, capture_output=True, text=True, timeout=60)
    except subprocess.TimeoutExpired:
        print(f"  TIMEOUT — skipping {name}")
        continue

    if os.path.exists(pdbqt_path):
        print(f"  PDBQT generated: {pdbqt_path}")
    else:
        print(f"  ERROR generating PDBQT: {result2.stderr}")

print("\n" + "=" * 55)
print("Ligand preparation complete.")

print("\nGenerated files:")
for f in sorted(os.listdir(ligands_dir)):
    size = os.path.getsize(os.path.join(ligands_dir, f))
    print(f"  {f} ({size} bytes)")