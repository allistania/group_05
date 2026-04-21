# -*- coding: utf-8 -*-
import subprocess
import re
import time
import os
import numpy as np
import matplotlib.pyplot as plt

# Slurm parameters (adjust to your cluster)
SBATCH_OPTS = {
    'partition': 'gpu',      # partition name
    'time': '00:11:00',          # time limit per run
    'cpus_per_task': 1,
    'nodelist': 'cluster4' 
}

EXE_PATH = "./shem"               # имя вашей MPI-программы
MESH_FILE = "fin2.msh"            # <-- добавьте имя файла сетки
MPI_MODULE = "openmpi/4.1.5"

def submit_job(n_procs, mesh_file, job_name=None):
    if job_name is None:
        job_name = f"mpi_p{n_procs}"
    
    out_file = f"{job_name}.out"
    err_file = f"{job_name}.err"
    
    script_content = f"""#!/bin/bash
#SBATCH -J {job_name}
#SBATCH --ntasks={n_procs}
#SBATCH --cpus-per-task={SBATCH_OPTS['cpus_per_task']}
#SBATCH --time={SBATCH_OPTS['time']}
#SBATCH --partition={SBATCH_OPTS['partition']}
#SBATCH --output={out_file}
#SBATCH --error={err_file}

module load {MPI_MODULE}
mpirun {EXE_PATH} {mesh_file}
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

def main():
    p_values = list(range(1, 10))
    results = []
    job_map = {}

    print("Submitting jobs...")
    for p in p_values:
        job_id, out_file, err_file = submit_job(p, MESH_FILE, job_name=f"dz1_p{p}")
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

    # Вывод таблицы
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

    # Построение графиков
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

if __name__ == "__main__":
    main()