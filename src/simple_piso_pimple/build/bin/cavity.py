import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

def read_output(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    data_lines = []
    header = None
    for line in lines:
        if line.startswith('#'):
            continue
        if header is None:
            header = line.strip().split(',')
            continue
        data_lines.append(line.strip())
    if header is None:
        header = ['x', 'y', 'rho', 'u', 'v', 'P', 'e_internal']
    df = pd.DataFrame([list(map(float, l.split(','))) for l in data_lines], columns=header)
    return df

df = read_output('output.csv')

x = df['x'].values
y = df['y'].values
u = df['u'].values
v = df['v'].values

nx = len(np.unique(x))
ny = len(np.unique(y))

X = x.reshape(ny, nx)
Y = y.reshape(ny, nx)
U = u.reshape(ny, nx)
V = v.reshape(ny, nx)

ix_center = np.argmin(np.abs(X[0,:] - 0.5))
iy_center = np.argmin(np.abs(Y[:,0] - 0.5))

ux_center = U[:, ix_center]
vy_center = V[iy_center, :]

y_line = Y[:, ix_center]
x_line = X[iy_center, :]

ghia_u_y = np.array([1.000, 0.977, 0.969, 0.953, 0.852, 0.734, 0.617, 0.500,
                     0.453, 0.281, 0.172, 0.000])
ghia_u_val = np.array([1.000, 0.841, 0.789, 0.687, 0.232, 0.003, -0.136, -0.206,
                       -0.211, -0.157, -0.102, 0.000])

ghia_v_x = np.array([1.000, 0.969, 0.961, 0.953, 0.852, 0.734, 0.617, 0.500,
                     0.453, 0.281, 0.172, 0.000])
ghia_v_val = np.array([0.000, -0.217, -0.332, -0.373, -0.232, -0.152, 0.058, 0.180,
                       0.172, 0.175, 0.162, 0.000])

interp_u = interp1d(y_line, ux_center, kind='linear', fill_value='extrapolate')
interp_v = interp1d(x_line, vy_center, kind='linear', fill_value='extrapolate')

u_calc = interp_u(ghia_u_y)
v_calc = interp_v(ghia_v_x)

u_err = np.abs((u_calc - ghia_u_val) / ghia_u_val) * 100
v_err = np.abs((v_calc - ghia_v_val) / ghia_v_val) * 100

with open('ghia_table.txt', 'w') as f:
    f.write('Comparison with Ghia (1982) Re=100\n')
    f.write('='*70 + '\n')
    f.write('Profile u(y) at x=0.5:\n')
    f.write(f"{'y':>8} {'u_calc':>10} {'u_ghia':>10} {'error %':>10}\n")
    for i in range(len(ghia_u_y)):
        f.write(f"{ghia_u_y[i]:8.4f} {u_calc[i]:10.4f} {ghia_u_val[i]:10.4f} {u_err[i]:9.2f}\n")
    f.write('\nProfile v(x) at y=0.5:\n')
    f.write(f"{'x':>8} {'v_calc':>10} {'v_ghia':>10} {'error %':>10}\n")
    for i in range(len(ghia_v_x)):
        f.write(f"{ghia_v_x[i]:8.4f} {v_calc[i]:10.4f} {ghia_v_val[i]:10.4f} {v_err[i]:9.2f}\n")

print('Table saved to ghia_table.txt')

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12,5))

ax1.plot(y_line, ux_center, 'b-', label='Computed')
ax1.plot(ghia_u_y, ghia_u_val, 'ro', label='Ghia (1982)', markersize=4)
ax1.set_xlabel('y')
ax1.set_ylabel('u')
ax1.set_title('u(y) at x=0.5, Re=100')
ax1.legend()
ax1.grid(True)

ax2.plot(x_line, vy_center, 'b-', label='Computed')
ax2.plot(ghia_v_x, ghia_v_val, 'ro', label='Ghia (1982)', markersize=4)
ax2.set_xlabel('x')
ax2.set_ylabel('v')
ax2.set_title('v(x) at y=0.5, Re=100')
ax2.legend()
ax2.grid(True)

plt.tight_layout()
plt.savefig('profiles.png', dpi=150)
plt.show()

dx = X[0,1] - X[0,0]
dy = Y[1,0] - Y[0,0]

vorticity = np.zeros_like(U)
for i in range(1, nx-1):
    for j in range(1, ny-1):
        dvdx = (V[j, i+1] - V[j, i-1]) / (2*dx)
        dudy = (U[j+1, i] - U[j-1, i]) / (2*dy)
        vorticity[j,i] = dvdx - dudy

psi = np.zeros_like(U)
psi[:,0] = 0.0
psi[0,:] = 0.0
psi[:,-1] = 0.0
psi[-1,:] = 0.0

max_iter = 10000
tol = 1e-6
for it in range(max_iter):
    psi_old = psi.copy()
    for i in range(1, nx-1):
        for j in range(1, ny-1):
            psi[j,i] = 0.25 * (psi[j, i-1] + psi[j, i+1] +
                               psi[j-1, i] + psi[j+1, i] +
                               dx*dy * vorticity[j,i])
    diff = np.max(np.abs(psi - psi_old))
    if diff < tol:
        break

plt.figure(figsize=(8,6))
levels = np.linspace(np.min(psi), np.max(psi), 30)
contour = plt.contour(X, Y, psi, levels=levels, cmap='viridis')
plt.colorbar(contour, label='Stream function')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Streamlines, Re=100')
plt.axis('equal')
plt.savefig('streamlines.png', dpi=150)
plt.show()

print(f'Converged in {it+1} iterations, residual={diff:.2e}')