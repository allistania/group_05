# -*- coding: utf-8 -*-
import subprocess
import re
import time
import os
import numpy as np
import matplotlib.pyplot as plt
import csv
import math

# Slurm parameters
SBATCH_OPTS = {
    'partition': 'gpu',
    'time': '00:11:00',
    'cpus_per_task': 1,
    'nodelist': 'cluster4'
}

EXE_PATH = "./shem"
MESH_FILE = "fin2.msh"
MPI_MODULE = "openmpi/4.1.5"

def submit_job(n_procs, mesh_file, prefix=None, job_name=None):
    if job_name is None:
        job_name = f"mpi_p{n_procs}"
    out_file = f"{job_name}.out"
    err_file = f"{job_name}.err"
    
    prefix_arg = f" {prefix}" if prefix else ""
    
    script_content = f"""#!/bin/bash
#SBATCH -J {job_name}
#SBATCH --ntasks={n_procs}
#SBATCH --cpus-per-task={SBATCH_OPTS['cpus_per_task']}
#SBATCH --time={SBATCH_OPTS['time']}
#SBATCH --partition={SBATCH_OPTS['partition']}
#SBATCH --output={out_file}
#SBATCH --error={err_file}

module load {MPI_MODULE}
mpirun {EXE_PATH} {mesh_file}{prefix_arg}
"""
    script_file = f"job_{job_name}.sh"
    with open(script_file, "w") as f:
        f.write(script_content)
    
    try:
        result = subprocess.run(["sbatch", script_file], capture_output=True, text=True, check=True)
        match = re.search(r"(\d+)", result.stdout)
        if match:
            job_id = int(match.group(1))
            print(f"Submitted job {job_id} for p={n_procs}")
            return job_id, out_file, err_file
        else:
            print(f"Failed to parse job id from: {result.stdout}")
            return None, None, None
    except subprocess.CalledProcessError as e:
        print(f"Error submitting job: {e.stderr}")
        return None, None, None
    finally:
        os.remove(script_file)

def wait_for_jobs(job_ids, check_interval=5):
    if not job_ids:
        return
    pending = set(str(jid) for jid in job_ids if jid is not None)
    while pending:
        time.sleep(check_interval)
        result = subprocess.run(["squeue", "-u", os.environ.get("USER", ""), "-h", "-o", "%i"],
                                capture_output=True, text=True)
        active = set(result.stdout.strip().split())
        pending = pending & active
        if pending:
            print(f"Waiting for jobs: {', '.join(pending)} ...")

def parse_time_from_file(out_file):
    try:
        with open(out_file, "r") as f:
            content = f.read()
        match = re.search(r"Total wall time:\s*([0-9.eE+-]+)", content)
        if match:
            return float(match.group(1))
        else:
            print(f"Total wall time not found in {out_file}")
            return None
    except FileNotFoundError:
        print(f"Output file {out_file} not found")
        return None

def parse_vtk(filename):
    rho = []
    vel = []
    try:
        with open(filename, "r") as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f"VTK file {filename} not found.")
        return None, None

    in_cell_data = False
    reading_density = False
    reading_velocity = False
    density_values = []
    vel_data = []

    for line in lines:
        stripped = line.strip()
        if stripped.startswith("CELL_DATA"):
            in_cell_data = True
            continue
        if not in_cell_data:
            continue

        if stripped.startswith("SCALARS") and "density" in stripped.lower():
            reading_density = True
            reading_velocity = False
            continue
        if stripped.startswith("VECTORS") and "velocity" in stripped.lower():
            reading_velocity = True
            reading_density = False
            continue
        if stripped.startswith("LOOKUP_TABLE"):
            continue

        if reading_density:
            parts = stripped.split()
            for p in parts:
                try:
                    density_values.append(float(p))
                except ValueError:
                    pass
            if len(density_values) > 0 and (stripped == "" or stripped.startswith("SCALARS") or stripped.startswith("VECTORS")):
                reading_density = False
                rho = density_values
                density_values = []

        if reading_velocity:
            parts = stripped.split()
            for i in range(0, len(parts), 3):
                if i+2 < len(parts):
                    try:
                        u = float(parts[i])
                        v = float(parts[i+1])
                        vel_data.append((u, v))
                    except ValueError:
                        pass
            if stripped == "" or stripped.startswith("SCALARS") or stripped.startswith("VECTORS"):
                reading_velocity = False
                vel = [math.sqrt(u*u + v*v) for (u,v) in vel_data]
                vel_data = []

    if reading_density and density_values:
        rho = density_values
    if reading_velocity and vel_data:
        vel = [math.sqrt(u*u + v*v) for (u,v) in vel_data]

    if not rho or not vel:
        print(f"Warning: Could not parse density or velocity from {filename}")
    return rho, vel

def compare_solutions(seq_file, par_file):
    rho_seq, u_seq = parse_vtk(seq_file)
    rho_par, u_par = parse_vtk(par_file)
    if rho_seq is None or rho_par is None or u_seq is None or u_par is None:
        return None, None
    if len(rho_seq) != len(rho_par) or len(u_seq) != len(u_par):
        print(f"Error: Mismatch in cell counts: seq={len(rho_seq)}, par={len(rho_par)}")
        return None, None

    max_rho_diff = max(abs(rho_par[i] - rho_seq[i]) for i in range(len(rho_seq)))
    max_u_diff   = max(abs(u_par[i] - u_seq[i]) for i in range(len(u_seq)))
    return max_rho_diff, max_u_diff

def main():
    p_values = list(range(1, 10))
    results = []
    job_map = {}

    print("Submitting jobs...")
    for p in p_values:
        job_id, out_file, err_file = submit_job(p, MESH_FILE, prefix=None, job_name=f"dz1_p{p}")
        if job_id:
            job_map[p] = (job_id, out_file, err_file)

    if not job_map:
        print("No jobs submitted. Exiting.")
        return

    print("\nWaiting for all jobs to complete...")
    wait_for_jobs([job_id for (job_id, _, _) in job_map.values()])

    print("\nCollecting results...")
    for p in p_values:
        if p in job_map:
            _, out_file, err_file = job_map[p]
            t = parse_time_from_file(out_file)
            results.append(t if t is not None else np.nan)
        else:
            results.append(np.nan)

    print("\n========== Performance Summary ==========")
    if not np.isnan(results[0]):
        serial_time = results[0]
        print(f"Serial time (p=1): {serial_time:.4f} s")
        print(f"{'p':>3} {'Time (s)':>10} {'Speedup':>10} {'Efficiency':>10}")
        for i, p in enumerate(p_values):
            t = results[i]
            if not np.isnan(t):
                speedup = serial_time / t
                efficiency = speedup / p
                print(f"{p:3d} {t:10.4f} {speedup:10.4f} {efficiency:10.4f}")
            else:
                print(f"{p:3d} {'N/A':>10} {'N/A':>10} {'N/A':>10}")
    else:
        print("No valid time for p=1, cannot compute speedup.")

    csv_filename = "performance.csv"
    with open(csv_filename, "w", newline="", encoding="utf-8") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["p", "Time (s)", "Speedup", "Efficiency"])
        if not np.isnan(results[0]):
            serial_time = results[0]
            for i, p in enumerate(p_values):
                t = results[i]
                if not np.isnan(t):
                    speedup = serial_time / t
                    efficiency = speedup / p
                    writer.writerow([p, f"{t:.6f}", f"{speedup:.6f}", f"{efficiency:.6f}"])
                else:
                    writer.writerow([p, "N/A", "N/A", "N/A"])
        else:
            for i, p in enumerate(p_values):
                writer.writerow([p, "N/A", "N/A", "N/A"])
    print(f"\nPerformance data saved to {csv_filename}")

    plt.figure(figsize=(12,5))
    plt.subplot(1,2,1)
    if not np.isnan(results[0]):
        speedup = [results[0]/t if not np.isnan(t) else np.nan for t in results]
        plt.plot(p_values, speedup, 'o-', label='Measured')
    plt.plot(p_values, p_values, 'k--', label='Ideal')
    plt.xlabel('Processes p')
    plt.ylabel('Speedup S')
    plt.title('Speedup')
    plt.legend()
    plt.grid(True)

    plt.subplot(1,2,2)
    if not np.isnan(results[0]):
        efficiency = [s/p if not np.isnan(s) else np.nan for s,p in zip(speedup, p_values)]
        plt.plot(p_values, efficiency, 's-', label='Measured')
    plt.plot(p_values, [1]*len(p_values), 'k--', label='Ideal')
    plt.xlabel('Processes p')
    plt.ylabel('Efficiency E')
    plt.title('Efficiency')
    plt.legend()
    plt.grid(True)

    plt.tight_layout()
    plt.savefig('performance.png')
    plt.show()

    print("\n--- Running verification (seq vs par) ---")
    seq_job_id, seq_out, seq_err = submit_job(1, MESH_FILE, prefix="seq", job_name="verify_seq")
    par_job_id, par_out, par_err = submit_job(8, MESH_FILE, prefix="par", job_name="verify_par")
    
    job_ids = []
    if seq_job_id: job_ids.append(seq_job_id)
    if par_job_id: job_ids.append(par_job_id)
    
    if job_ids:
        print("Waiting for verification jobs to complete...")
        wait_for_jobs(job_ids)
    
    seq_vtk = "seq_final.vtk"
    par_vtk = "par_final.vtk"
    max_rho, max_u = compare_solutions(seq_vtk, par_vtk)
    
    if max_rho is not None and max_u is not None:
        with open("verification.txt", "w") as vf:
            vf.write(f"max|rho_par-rho_seq| = {max_rho:.1e}\n")
            vf.write(f"max|u_par-u_seq| = {max_u:.1e}\n")
        print(f"Verification results saved to verification.txt")
        print(f"max|rho_par-rho_seq| = {max_rho:.1e}")
        print(f"max|u_par-u_seq| = {max_u:.1e}")
    else:
        print("Verification failed: could not compare VTK files.")

if __name__ == "__main__":
    main()