# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd
from matplotlib.ticker import LogLocator

# Create directory for plots
os.makedirs('plots', exist_ok=True)

# Initialize variables to avoid reference errors
t_start1, xs_start1, rhos_start1, Ps_start1, us_start1, Is_start1 = [], [], [], [], [], []

# Try to read first dataset (exact solution)
try:
    with open("data/tochnoeresh.csv") as f:
        labels = f.readline()
        data1 = f.readlines()
    
    data1 = list(map(lambda x: x.split(','), data1))
    
    t_start1 = list(map(lambda x: float(x[0]), data1))
    xs_start1 = list(map(lambda x: float(x[1]), data1))
    rhos_start1 = list(map(lambda x: float(x[2]), data1))
    Ps_start1 = list(map(lambda x: float(x[3]), data1))
    us_start1 = list(map(lambda x: float(x[4]), data1))
    Is_start1 = list(map(lambda x: float(x[5]), data1))
    print("Successfully loaded tochnoeresh.csv")
except FileNotFoundError:
    print("Warning: File tochnoeresh.csv not found")
except Exception as e:
    print(f"Error reading tochnoeresh.csv: {e}")

# Read other datasets including Custom
datasets = [
    "data/hislennoeresh_tochnoe.csv",
    "data/hislennoeresh_hll.csv", 
    "data/hislennoeresh_hllc.csv",
    "data/hislennoeresh_roe.csv",
    "data/hislennoeresh_osher.csv",
    "data/hislennoereshKabare.csv",
    "data/hislennoeresh_rusanov.csv"  # Добавлен Custom файл
]

data_list = []
for dataset in datasets:
    try:
        with open(dataset) as f:
            labels = f.readline()
            data = f.readlines()
        data = list(map(lambda x: x.split(','), data))
        data_list.append(data)
        print(f"Successfully loaded {dataset}")
    except FileNotFoundError:
        print(f"Warning: File {dataset} not found")
        data_list.append([])
    except Exception as e:
        print(f"Error reading {dataset}: {e}")
        data_list.append([])

# Extract data for each dataset
t_starts = []
xs_starts = []
rhos_starts = []
Ps_starts = []
us_starts = []
Is_starts = []

for data in data_list:
    if data:  # Check that data exists
        try:
            t_starts.append(list(map(lambda x: float(x[0]), data)))
            xs_starts.append(list(map(lambda x: float(x[1]), data)))
            rhos_starts.append(list(map(lambda x: float(x[2]), data)))
            Ps_starts.append(list(map(lambda x: float(x[3]), data)))
            us_starts.append(list(map(lambda x: float(x[4]), data)))
            Is_starts.append(list(map(lambda x: float(x[5]), data)))
        except Exception as e:
            print(f"Error processing data: {e}")
            # Add empty lists if processing fails
            t_starts.append([])
            xs_starts.append([])
            rhos_starts.append([])
            Ps_starts.append([])
            us_starts.append([])
            Is_starts.append([])
    else:
        # Add empty lists if file not found
        t_starts.append([])
        xs_starts.append([])
        rhos_starts.append([])
        Ps_starts.append([])
        us_starts.append([])
        Is_starts.append([])

# Check if we have any data to plot
has_exact_solution = len(xs_start1) > 0
has_numerical_solutions = any(len(xs) > 0 for xs in xs_starts)

if not has_exact_solution and not has_numerical_solutions:
    print("Error: No data files found to plot!")
    exit(1)

# Determine used_time
used_time = None
if has_exact_solution:
    used_time = t_start1[0]
elif has_numerical_solutions:
    # Find first non-empty numerical solution for time reference
    for i in range(len(t_starts)):
        if t_starts[i]:
            used_time = t_starts[i][0]
            break

if used_time is None:
    print("Error: Could not determine time for plot!")
    exit(1)

# Verify time matches (only if we have both exact and numerical solutions)
if has_exact_solution and has_numerical_solutions:
    time_mismatch = False
    for i in range(len(t_starts)):
        if t_starts[i] and abs(t_start1[0] - t_starts[i][0]) > 1e-10:
            time_mismatch = True
            print(f"Warning: Time mismatch! t1={t_start1[0]}, t{i+1}={t_starts[i][0]}")

# Colors and labels for different schemes (добавлен Custom)
colors = ['#086522', '#FF4500', '#00FFFF', '#FF00FF', '#0000FF', '#00FF00', '#FFA500']
labels = [
    "Exact solution", 
    "Riehman", 
    "HLL", 
    "HLLC",
    "ROE",
    "OSHER",
    "Kabare scheme",
    "RUSANOV"
]

markers = ['', 'o', 'o', 'o', 'o', 'o', 'o']
markersizes = [0, 2, 2, 2, 2, 2, 2]

# Create the main combined plot with all 4 variables
fig, axs = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle(f"Solve GD equations at t={used_time:.3f} s.", fontsize=14)

# Pressure plot
axs[0,0].grid(True, linestyle='--', alpha=0.7)
axs[0,0].set_xlabel("x")
axs[0,0].set_ylabel("P")
if has_exact_solution:
    axs[0,0].plot(xs_start1, Ps_start1, color=colors[0], label=labels[0], linewidth=2)
for i in range(len(datasets)):
    if xs_starts[i] and Ps_starts[i]:
        marker = markers[i+1] if i+1 < len(markers) else 'o'
        axs[0,0].plot(xs_starts[i], Ps_starts[i], '-'+marker, color=colors[i+1], 
                      markersize=markersizes[i+1], label=labels[i+1], linewidth=1)
if has_exact_solution or has_numerical_solutions:
    axs[0,0].legend()

# Velocity plot
axs[0,1].grid(True, linestyle='--', alpha=0.7)
axs[0,1].set_xlabel("x")
axs[0,1].set_ylabel("U")
if has_exact_solution:
    axs[0,1].plot(xs_start1, us_start1, color=colors[0], label=labels[0], linewidth=2)
for i in range(len(datasets)):
    if xs_starts[i] and us_starts[i]:
        marker = markers[i+1] if i+1 < len(markers) else 'o'
        axs[0,1].plot(xs_starts[i], us_starts[i], '-'+marker, color=colors[i+1], 
                      markersize=markersizes[i+1], label=labels[i+1], linewidth=1)
if has_exact_solution or has_numerical_solutions:
    axs[0,1].legend()

# Density plot
axs[1,0].grid(True, linestyle='--', alpha=0.7)
axs[1,0].set_xlabel("x")
axs[1,0].set_ylabel("rho")
if has_exact_solution:
    axs[1,0].plot(xs_start1, rhos_start1, color=colors[0], label=labels[0], linewidth=2)
for i in range(len(datasets)):
    if xs_starts[i] and rhos_starts[i]:
        marker = markers[i+1] if i+1 < len(markers) else 'o'
        axs[1,0].plot(xs_starts[i], rhos_starts[i], '-'+marker, color=colors[i+1], 
                      markersize=markersizes[i+1], label=labels[i+1], linewidth=1)
if has_exact_solution or has_numerical_solutions:
    axs[1,0].legend()

# Internal energy plot
axs[1,1].grid(True, linestyle='--', alpha=0.7)
axs[1,1].set_xlabel("x")
axs[1,1].set_ylabel("I")
if has_exact_solution:
    axs[1,1].plot(xs_start1, Is_start1, color=colors[0], label=labels[0], linewidth=2)
for i in range(len(datasets)):
    if xs_starts[i] and Is_starts[i]:
        marker = markers[i+1] if i+1 < len(markers) else 'o'
        axs[1,1].plot(xs_starts[i], Is_starts[i], '-'+marker, color=colors[i+1], 
                      markersize=markersizes[i+1], label=labels[i+1], linewidth=1)
if has_exact_solution or has_numerical_solutions:
    axs[1,1].legend()

plt.tight_layout()
fig.savefig('plots/combined_plot.png', dpi=300, bbox_inches='tight')
plt.close(fig)
print("Main plot saved as 'plots/combined_plot.png'")

# Create error plot if error file exists
try:
    # Read error data - assuming it's in the format you provided
    error_data = pd.read_csv("data/errL1.csv")
    print("Successfully loaded error data from data/errL1.csv")
    
    # Create error plot with 4 subplots
    fig, axs = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle("Error Convergence Analysis", fontsize=14)
    
    # Define colors and markers for error plots
    error_colors = ['blue', 'red', 'green', 'orange']
    error_markers = ['o', 's', '^', 'd']
    
    # Density error plot
    if 'Rho' in error_data.columns:
        axs[0,0].plot(error_data['X'], error_data['Rho'], 
                     color=error_colors[0], marker=error_markers[0],
                     linewidth=2, markersize=6, label='Density (Rho)')
        axs[0,0].set_xlabel("Grid spacing X")
        axs[0,0].set_ylabel("Rho Error")
        axs[0,0].set_title("Density Error Convergence")
        axs[0,0].set_xscale('log')
        axs[0,0].set_yscale('log')
        
        # Более частая сетка для логарифмических шкал
        axs[0,0].grid(True, which='both', linestyle=':', alpha=0.3)
        axs[0,0].grid(True, which='major', linestyle='-', alpha=0.5)
        
        # Устанавливаем более частые деления для сетки
        axs[0,0].xaxis.set_major_locator(LogLocator(numticks=15))
        axs[0,0].yaxis.set_major_locator(LogLocator(numticks=15))
        
        axs[0,0].legend()
    
    # Velocity error plot
    if 'U' in error_data.columns:
        axs[0,1].plot(error_data['X'], error_data['U'], 
                     color=error_colors[1], marker=error_markers[1],
                     linewidth=2, markersize=6, label='Velocity (U)')
        axs[0,1].set_xlabel("Grid spacing X")
        axs[0,1].set_ylabel("U Error")
        axs[0,1].set_title("Velocity Error Convergence")
        axs[0,1].set_xscale('log')
        axs[0,1].set_yscale('log')
        
        # Более частая сетка для логарифмических шкал
        axs[0,1].grid(True, which='both', linestyle=':', alpha=0.3)
        axs[0,1].grid(True, which='major', linestyle='-', alpha=0.5)
        
        # Устанавливаем более частые деления для сетки
        axs[0,1].xaxis.set_major_locator(LogLocator(numticks=15))
        axs[0,1].yaxis.set_major_locator(LogLocator(numticks=15))
        
        axs[0,1].legend()
    
    # Pressure error plot
    if 'P' in error_data.columns:
        axs[1,0].plot(error_data['X'], error_data['P'], 
                     color=error_colors[2], marker=error_markers[2],
                     linewidth=2, markersize=6, label='Pressure (P)')
        axs[1,0].set_xlabel("Grid spacing X")
        axs[1,0].set_ylabel("P Error")
        axs[1,0].set_title("Pressure Error Convergence")
        axs[1,0].set_xscale('log')
        axs[1,0].set_yscale('log')
        
        # Более частая сетка для логарифмических шкал
        axs[1,0].grid(True, which='both', linestyle=':', alpha=0.3)
        axs[1,0].grid(True, which='major', linestyle='-', alpha=0.5)
        
        # Устанавливаем более частые деления для сетки
        axs[1,0].xaxis.set_major_locator(LogLocator(numticks=15))
        axs[1,0].yaxis.set_major_locator(LogLocator(numticks=15))
        
        axs[1,0].legend()
    
    # Internal energy error plot
    if 'I' in error_data.columns:
        axs[1,1].plot(error_data['X'], error_data['I'], 
                     color=error_colors[3], marker=error_markers[3],
                     linewidth=2, markersize=6, label='Internal Energy (I)')
        axs[1,1].set_xlabel("Grid spacing X")
        axs[1,1].set_ylabel("I Error")
        axs[1,1].set_title("Internal Energy Error Convergence")
        axs[1,1].set_xscale('log')
        axs[1,1].set_yscale('log')
        
        # Более частая сетка для логарифмических шкал
        axs[1,1].grid(True, which='both', linestyle=':', alpha=0.3)
        axs[1,1].grid(True, which='major', linestyle='-', alpha=0.5)
        
        # Устанавливаем более частые деления для сетки
        axs[1,1].xaxis.set_major_locator(LogLocator(numticks=15))
        axs[1,1].yaxis.set_major_locator(LogLocator(numticks=15))
        
        axs[1,1].legend()
    
    plt.tight_layout()
    
    # Save error plot
    fig.savefig('plots/error_convergence.png', dpi=300, bbox_inches='tight')
    plt.close(fig)
    print("Error plot saved as 'plots/error_convergence.png'")
    
except FileNotFoundError:
    print("Warning: Error file data/errL1.csv not found")
except Exception as e:
    print(f"Error processing error file: {e}")

print("Script completed successfully!")