import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# -------------------------------------------------------------------
# 1. Reading data from file (skip comments AND header)
# -------------------------------------------------------------------
filename = "output.csv"

data = []
with open(filename, 'r') as f:
    # Skip all lines starting with '#'
    for line in f:
        if line.startswith('#'):
            continue
        else:
            # First non-comment line is the header: "x,y,rho,u,v,P,e_internal"
            # We skip it by breaking and then reading the rest
            break
    # Now read the actual data
    for line in f:
        line = line.strip()
        if not line:
            continue
        parts = line.split(',')
        data.append([float(x) for x in parts])

data = np.array(data)
gamma = 1.4
x = data[:, 0]
y = data[:, 1]
rho = data[:, 2]
u = data[:, 3]
v = data[:, 4]
P = data[:, 5]
e_int = P/(rho*(gamma-1))
#e_int = data[:, 6]                     # internal energy from file

# -------------------------------------------------------------------
# 2. Grid dimensions and reshaping
# -------------------------------------------------------------------
nx = len(np.unique(x))
ny = len(np.unique(y))

X = x.reshape(ny, nx)
Y = y.reshape(ny, nx)
Rho = rho.reshape(ny, nx)
U = u.reshape(ny, nx)
V = v.reshape(ny, nx)
P = P.reshape(ny, nx)
E_int = P/(Rho*(gamma-1))
#E_int = e_int.reshape(ny, nx)

# Velocity magnitude
speed = np.sqrt(U**2 + V**2)

# -------------------------------------------------------------------
# 3. Plotting – 2?3 panel with all variables
# -------------------------------------------------------------------
plt.style.use('seaborn-v0_8-whitegrid')
fig, axes = plt.subplots(2, 3, figsize=(16, 11))
fig.suptitle('2D Euler solver – flow fields', fontsize=16, y=0.98)

# ---- Density ----
ax = axes[0, 0]
im = ax.pcolormesh(X, Y, Rho, shading='auto', cmap='viridis')
ax.set_title(r'Density $\rho$', fontsize=12)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_aspect('equal')
plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

# ---- Pressure ----
ax = axes[0, 1]
im = ax.pcolormesh(X, Y, P, shading='auto', cmap='plasma')
ax.set_title(r'Pressure $P$', fontsize=12)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_aspect('equal')
plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

# ---- Velocity u ----
ax = axes[0, 2]
im = ax.pcolormesh(X, Y, U, shading='auto', cmap='coolwarm')
ax.set_title(r'Velocity $u$', fontsize=12)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_aspect('equal')
plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

# ---- Velocity v ----
ax = axes[1, 0]
im = ax.pcolormesh(X, Y, V, shading='auto', cmap='coolwarm')
ax.set_title(r'Velocity $v$', fontsize=12)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_aspect('equal')
plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

# ---- Internal energy ----
ax = axes[1, 1]
im = ax.pcolormesh(X, Y, E_int, shading='auto', cmap='inferno')
ax.set_title(r'Internal energy $e_{int}$', fontsize=12)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_aspect('equal')
plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

# ---- Velocity magnitude ----
ax = axes[1, 2]
im = ax.pcolormesh(X, Y, speed, shading='auto', cmap='magma')
ax.set_title(r'Velocity magnitude $\sqrt{u^2+v^2}$', fontsize=12)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_aspect('equal')
plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

plt.tight_layout()
plt.savefig('plots/godunov2D_fields.png', dpi=150, bbox_inches='tight')
plt.show()

# -------------------------------------------------------------------
# 4. 3D density surface
# -------------------------------------------------------------------
fig3d = plt.figure(figsize=(10, 8))
ax3d = fig3d.add_subplot(111, projection='3d')
surf = ax3d.plot_surface(X, Y, Rho, cmap='viridis', edgecolor='none', alpha=0.9)
ax3d.set_title(r'3D surface of density $\rho(x,y)$', fontsize=14)
ax3d.set_xlabel('x')
ax3d.set_ylabel('y')
ax3d.set_zlabel(r'$\rho$')
fig3d.colorbar(surf, ax=ax3d, shrink=0.5, aspect=10)

plt.tight_layout()
plt.savefig('plots/godunov2D_3d_rho.png', dpi=150)
plt.show()

# -------------------------------------------------------------------
# 5. Cross sections at x=0.5 and y=0.5 for all variables
# -------------------------------------------------------------------
variables = {
    'rho': r'Density $\rho$',
    'u': r'Velocity $u$',
    'v': r'Velocity $v$',
    'P': r'Pressure $P$',
    'e_int': r'Internal energy $e_{int}$',
    'speed': r'Velocity magnitude'
}
data_vars = [rho, u, v, P, e_int, speed]
var_labels = list(variables.values())

unique_x = np.unique(x)
unique_y = np.unique(y)
ix_center = np.argmin(np.abs(unique_x - 0.5))
iy_center = np.argmin(np.abs(unique_y - 0.5))

fig_profiles, axes = plt.subplots(len(var_labels), 2, figsize=(14, 3*len(var_labels)))
fig_profiles.suptitle('Cross sections at x = 0.5 (left) and y = 0.5 (right)', fontsize=16, y=0.98)

for idx, (var_data, label) in enumerate(zip(data_vars, var_labels)):
    # reshape to grid
    Var = var_data.reshape(ny, nx)
    
    # profile at x=0.5 (constant x) -> along y
    axes[idx, 0].plot(unique_y, Var[:, ix_center], 'b-', linewidth=1.5, markersize=3)
    axes[idx, 0].set_xlabel('y')
    axes[idx, 0].set_ylabel(label)
    axes[idx, 0].grid(True, linestyle=':', alpha=0.7)
    
    # profile at y=0.5 (constant y) -> along x
    axes[idx, 1].plot(unique_x, Var[iy_center, :], 'r-', linewidth=1.5, markersize=3)
    axes[idx, 1].set_xlabel('x')
    axes[idx, 1].set_ylabel(label)
    axes[idx, 1].grid(True, linestyle=':', alpha=0.7)

plt.tight_layout()
plt.savefig('plots/godunov2D_all_profiles.png', dpi=150, bbox_inches='tight')
plt.show()