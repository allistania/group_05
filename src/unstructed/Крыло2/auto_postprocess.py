import os
import glob
import subprocess
import re

def auto_postprocess():
    results_dir = "results"
    pattern = os.path.join(results_dir, "*_cp_distribution.dat")
    for cp_file in glob.glob(pattern):
        base = os.path.basename(cp_file).replace("_cp_distribution.dat", "")
        quality_csv = os.path.join(results_dir, "mesh_quality", f"{base}_quality_raw.csv")
        if os.path.exists(quality_csv):
            print(f"Processing {base}...")
            subprocess.run(["python3", "Cp.py", cp_file, f"{base}_cp.png"])
            subprocess.run(["python3", "mesh_quality.py", quality_csv, base])
        else:
            print(f"Skipping {base}: quality CSV not found")
    print("All done.")

if __name__ == "__main__":
    auto_postprocess()