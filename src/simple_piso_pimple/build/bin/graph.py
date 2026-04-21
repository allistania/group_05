import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import griddata

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
speed_mag = np.sqrt(u**2 + v**2)

# -------------------------------------------------------------------
# 2. Original grid dimensions and reshaping
# -------------------------------------------------------------------
nx_orig = len(np.unique(x))
ny_orig = len(np.unique(y))

X_orig = x.reshape(ny_orig, nx_orig)
Y_orig = y.reshape(ny_orig, nx_orig)
Rho = rho.reshape(ny_orig, nx_orig)
U_orig = u.reshape(ny_orig, nx_orig)
V_orig = v.reshape(ny_orig, nx_orig)
P_orig = P.reshape(ny_orig, nx_orig)
E_int = P_orig/(Rho*(gamma-1))
speed = np.sqrt(U_orig**2 + V_orig**2)

# -------------------------------------------------------------------
# 3. Create uniform grid for streamlines
# -------------------------------------------------------------------
# Determine number of points for uniform grid (can be adjusted)
nx_uniform = 100
ny_uniform = 100

x_min, x_max = x.min(), x.max()
y_min, y_max = y.min(), y.max()

X_uniform = np.linspace(x_min, x_max, nx_uniform)
Y_uniform = np.linspace(y_min, y_max, ny_uniform)
X_grid, Y_grid = np.meshgrid(X_uniform, Y_uniform)

# Interpolate velocity components onto uniform grid
points = np.column_stack((x, y))
U_uniform = griddata(points, u, (X_grid, Y_grid), method='linear')
V_uniform = griddata(points, v, (X_grid, Y_grid), method='linear')

# Optional: handle NaN values (e.g., fill with 0 or nearest neighbor)
U_uniform = np.nan_to_num(U_uniform, nan=0.0)
V_uniform = np.nan_to_num(V_uniform, nan=0.0)

# -------------------------------------------------------------------
# 4. Plotting with streamlines on uniform grid
# -------------------------------------------------------------------
plt.style.use('seaborn-v0_8-whitegrid')
fig = plt.figure(figsize=(14, 10))

# Helper to add streamlines to any axis (uses uniform grid)
def add_streamlines(ax, X_grid, Y_grid, U_uniform, V_uniform, density=1.5):
    ax.streamplot(X_grid, Y_grid, U_uniform, V_uniform,
                  color='white', linewidth=0.8, density=density)

# ---- Density ----
ax1 = fig.add_subplot(2, 3, 1)
im1 = ax1.pcolormesh(X_orig, Y_orig, Rho, shading='auto', cmap='viridis')
add_streamlines(ax1, X_grid, Y_grid, U_uniform, V_uniform)
ax1.set_title(r'Density $\rho$ + streamlines')
ax1.set_xlabel('x')
ax1.set_ylabel('y')
ax1.set_aspect('equal')
plt.colorbar(im1, ax=ax1)

# ---- Pressure ----
ax2 = fig.add_subplot(2, 3, 2)
im2 = ax2.pcolormesh(X_orig, Y_orig, P_orig, shading='auto', cmap='viridis')
add_streamlines(ax2, X_grid, Y_grid, U_uniform, V_uniform)
ax2.set_title(r'Pressure $P$ + streamlines')
ax2.set_xlabel('x')
ax2.set_ylabel('y')
ax2.set_aspect('equal')
plt.colorbar(im2, ax=ax2)

# ---- Velocity u ----
ax3 = fig.add_subplot(2, 3, 3)
im3 = ax3.pcolormesh(X_orig, Y_orig, U_orig, shading='auto', cmap='viridis')
add_streamlines(ax3, X_grid, Y_grid, U_uniform, V_uniform)
ax3.set_title(r'Velocity $u$ + streamlines')
ax3.set_xlabel('x')
ax3.set_ylabel('y')
ax3.set_aspect('equal')
plt.colorbar(im3, ax=ax3)

# ---- Velocity v ----
ax4 = fig.add_subplot(2, 3, 4)
im4 = ax4.pcolormesh(X_orig, Y_orig, V_orig, shading='auto', cmap='viridis')
add_streamlines(ax4, X_grid, Y_grid, U_uniform, V_uniform)
ax4.set_title(r'Velocity $v$ + streamlines')
ax4.set_xlabel('x')
ax4.set_ylabel('y')
ax4.set_aspect('equal')
plt.colorbar(im4, ax=ax4)

# ---- Internal energy ----
ax5 = fig.add_subplot(2, 3, 5)
im5 = ax5.pcolormesh(X_orig, Y_orig, E_int, shading='auto', cmap='viridis')
add_streamlines(ax5, X_grid, Y_grid, U_uniform, V_uniform)
ax5.set_title(r'Internal energy $e_{int}$ + streamlines')
ax5.set_xlabel('x')
ax5.set_ylabel('y')
ax5.set_aspect('equal')
plt.colorbar(im5, ax=ax5)

# ---- Velocity magnitude ----
ax6 = fig.add_subplot(2, 3, 6)
im6 = ax6.pcolormesh(X_orig, Y_orig, speed, shading='auto', cmap='viridis')
add_streamlines(ax6, X_grid, Y_grid, U_uniform, V_uniform)
ax6.set_title(r'Velocity magnitude $\sqrt{u^2+v^2}$ + streamlines', fontsize=12)
ax6.set_xlabel('x')
ax6.set_ylabel('y')
ax6.set_aspect('equal')
plt.colorbar(im6, ax=ax6)

plt.tight_layout()
plt.savefig('plots/godunov2D_fields_with_streamlines.png', dpi=150)
plt.show()

# -------------------------------------------------------------------
# 5. 3D density surface (unchanged)
# -------------------------------------------------------------------
fig3d = plt.figure(figsize=(10, 8))
ax3d = fig3d.add_subplot(111, projection='3d')
surf = ax3d.plot_surface(X_orig, Y_orig, Rho, cmap='viridis', edgecolor='none', alpha=0.9)
ax3d.set_title(r'3D surface of density $\rho(x,y)$')
ax3d.set_xlabel('x')
ax3d.set_ylabel('y')
ax3d.set_zlabel(r'$\rho$')
fig3d.colorbar(surf, ax=ax3d, shrink=0.5, aspect=10)

plt.tight_layout()
plt.savefig('plots/godunov2D_3d_rho.png', dpi=150)
plt.show()

# -------------------------------------------------------------------
# 6. Cross sections (unchanged)
# -------------------------------------------------------------------
variables = {
    'rho': r'Density $\rho$',
    'u': r'Velocity $u$',
    'v': r'Velocity $v$',
    'P': r'Pressure $P$',
    'e_int': r'Internal energy $e_{int}$',
    'speed': r'Velocity magnitude $\sqrt{u^2+v^2}$'
}
data_vars = [rho, u, v, P, e_int, speed_mag]
var_labels = list(variables.values())

unique_x = np.unique(x)
unique_y = np.unique(y)
ix_center = np.argmin(np.abs(unique_x - 0.5))
iy_center = np.argmin(np.abs(unique_y - 0.5))

fig_profiles, axes = plt.subplots(len(var_labels), 2, figsize=(12, 3*len(var_labels)))
fig_profiles.suptitle('Cross sections at x=0.5 (left) and y=0.5 (right)', fontsize=14)

for idx, (var_data, label) in enumerate(zip(data_vars, var_labels)):
    Var = var_data.reshape(ny_orig, nx_orig)
    
    axes[idx, 0].plot(unique_y, Var[:, ix_center], 'b-o', markersize=3)
    axes[idx, 0].set_xlabel('y')
    axes[idx, 0].set_ylabel(label)
    axes[idx, 0].grid(True)
    
    axes[idx, 1].plot(unique_x, Var[iy_center, :], 'r-o', markersize=3)
    axes[idx, 1].set_xlabel('x')
    axes[idx, 1].set_ylabel(label)
    axes[idx, 1].grid(True)

plt.tight_layout()
plt.savefig('plots/godunov2D_all_profiles.png', dpi=150)
plt.show()