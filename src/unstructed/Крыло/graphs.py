import os
import glob
import numpy as np
import matplotlib.pyplot as plt

# Physical Constants (matching your Naca.cpp)
GAMMA = 1.4
P_INF = 1.0
RHO_INF = 1.0
U_INF = 15.0  # Speed of sound is approx 1.18 in these units if P=1, Rho=1
CHORD = 1.0

def parse_vtk_full(file_path):
    """Parses coordinates and cell data (pressure, density, momentum) from VTK."""
    points = []
    data = {}
    
    with open(file_path, 'r') as f:
        lines = f.readlines()
        
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith("POINTS"):
            num_pts = int(line.split()[1])
            for j in range(1, num_pts + 1):
                points.append([float(val) for val in lines[i+j].split()])
            i += num_pts
        
        # Look for Cell Data or Point Data
        if "SCALARS" in line or "VECTORS" in line:
            var_name = line.split()[1]
            i += 2 # Skip lookup table
            values = []
            while len(values) < len(points):
                values.extend([float(val) for val in lines[i].split()])
                i += 1
            data[var_name] = np.array(values)
            continue
        i += 1
    
    return np.array(points), data

def plot_aerodynamics():
    print("=== Aerodynamic Analysis: Mach Map and Cp Plot ===")
    
    vtk_files = glob.glob("results/simulation/*.vtk")
    if not vtk_files:
        print("Error: No VTK files found.")
        return
    
    latest_vtk = sorted(vtk_files)[-1]
    pts, data = parse_vtk_full(latest_vtk)
    
    x, y = pts[:, 0], pts[:, 1]
    rho = data.get('density')
    p = data.get('pressure')
    # Assuming momentum components are stored as separate scalars or a vector
    # In your Naca.cpp, check if it's 'momentum_x' or 'velocity'
    m_x = data.get('momentum_x', np.zeros_like(rho))
    m_y = data.get('momentum_y', np.zeros_like(rho))

    # 1. Calculate Mach Number
    # velocity = momentum / density
    v_mag = np.sqrt(m_x**2 + m_y**2) / rho
    # sound speed a = sqrt(gamma * P / rho)
    a = np.sqrt(GAMMA * p / rho)
    mach = v_mag / a

    # 2. Calculate Cp
    q_inf = 0.5 * RHO_INF * (U_INF ** 2)
    cp = (p - P_INF) / q_inf

    # --- PLOT 1: MACH NUMBER MAP ---
    plt.figure(figsize=(12, 5))
    sc = plt.scatter(x, y, c=mach, cmap='jet', s=1)
    plt.colorbar(sc, label='Mach Number')
    plt.title(f'Mach Number Field (Final Step)')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.axis('equal')
    plt.xlim(-0.5, 1.5)
    plt.ylim(-0.5, 0.5)
    plt.savefig('results/mach_map.png', dpi=300)
    print("Saved: results/mach_map.png")

    # --- PLOT 2: Cp DISTRIBUTION ---
    # Filter points near the airfoil surface
    mask = (x >= 0.0) & (x <= 1.0) & (np.abs(y) <= 0.1)
    x_s, y_s, cp_s = x[mask], y[mask], cp[mask]
    
    plt.figure(figsize=(10, 6))
    up = y_s >= 0
    lo = y_s < 0
    
    idx_u = np.argsort(x_s[up])
    plt.plot(x_s[up][idx_u], cp_s[up][idx_u], 'b-', label='Upper Surface')
    
    idx_l = np.argsort(x_s[lo])
    plt.plot(x_s[lo][idx_l], cp_s[lo][idx_l], 'r--', label='Lower Surface')

    plt.gca().invert_yaxis()
    plt.title('Pressure Coefficient $C_p$ Distribution')
    plt.xlabel('$x/c$')
    plt.ylabel('$C_p$')
    plt.grid(True)
    plt.legend()
    plt.savefig('results/cp_distribution.png')
    print("Saved: results/cp_distribution.png")

    # 3. Print Coefficients if file exists
    coeff_file = glob.glob("results/*_coefficients.dat")
    if coeff_file:
        print("\nSummary from .dat file:")
        with open(coeff_file[0], 'r') as f:
            print(f.read())

if __name__ == "__main__":
    if not os.path.exists('results'):
        os.makedirs('results')
    plot_aerodynamics()