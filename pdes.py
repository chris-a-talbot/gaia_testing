import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
from noise import pnoise2
from tqdm import tqdm
import os

os.makedirs('temperate_forest_simulation', exist_ok=True)

# Parameters
grid_size = (128, 128)
dx = 1.0
dt = 0.1
total_time = 100
D_A = 0.05
v = (0.0, 0.2)
s_max = 0.15
K0 = 200
urban_frac = 0.25
min_population = 1e-3

# Initialize environment
np.random.seed(42)
x, y = np.mgrid[0:grid_size[0], 0:grid_size[1]]
y_norm = y / grid_size[1]

# Temperature field
T0 = 10 + 15 * (1 - y_norm)


def generate_conductance(seed):
    scale = 0.08
    conductance = np.zeros(grid_size)
    for i in range(grid_size[0]):
        for j in range(grid_size[1]):
            conductance[i, j] = pnoise2(i * scale, j * scale, octaves=6,
                                        persistence=0.6, base=seed)
    conductance = (conductance - conductance.min()) / (conductance.max() - conductance.min())
    conductance = 0.3 + 0.5 * conductance

    # FIXED: Proper urban/corridor assignment
    urban_mask = np.random.rand(*grid_size) < urban_frac
    corridor_prob = np.random.rand(*grid_size)
    conductance[urban_mask] = np.where(corridor_prob[urban_mask] < 0.1, 0.4, 0.05)

    return conductance


C = generate_conductance(42)

# Initial conditions
A = np.zeros(grid_size)
A[y_norm > 0.7] = 0.8
A = gaussian_filter(A, sigma=3)
A = np.clip(A, 0, 1)

N = np.zeros(grid_size) + K0 * (0.6 + 0.4 * np.random.rand(*grid_size))
N = np.clip(N, min_population, None)
N[C < 0.1] = min_population


def laplacian(f, C):
    flux = (
                   C * (np.roll(f, 1, axis=0) - f) +
                   C * (np.roll(f, -1, axis=0) - f) +
                   C * (np.roll(f, 1, axis=1) - f) +
                   C * (np.roll(f, -1, axis=1) - f)
           ) / dx ** 2
    return flux


def advection(f, vy):
    return vy * (np.roll(f, 1, axis=0) - f) / dx


# Main simulation loop
frame = 0
for t in tqdm(range(int(total_time / dt))):
    T = T0 + 0.015 * t * dt
    K = K0 * np.exp(-0.5 * ((T - 22) / 7) ** 2) * (C > 0.1)
    K = np.clip(K, min_population, None)
    s = s_max * np.exp(-0.5 * ((T - 25) / 10) ** 2)

    # Allele frequency update
    dA_diff = laplacian(A, C)
    dA_adv = advection(A, v[1])
    dA_sel = s * A * (1 - A)
    A_new = A + dt * (dA_diff + dA_adv + dA_sel) + np.random.normal(0, 0.005, grid_size)
    A = np.clip(A_new, 0, 1)

    # Population density update
    dN_growth = 0.1 * N * np.clip(1 - N / K, -0.5, 1)
    dN_diff = laplacian(N, C)
    dN_adv = advection(N, v[1])
    N_new = N + dt * (dN_growth + dN_diff + dN_adv)
    N = np.clip(N_new, min_population, 1.2 * K)
    N = np.where(N < min_population, min_population, N)

    # Visualization every 5 years
    if t % 50 == 0:
        fig, ax = plt.subplots(1, 3, figsize=(20, 6))

        # Allele frequency
        im0 = ax[0].imshow(A.T, cmap='RdYlBu', vmin=0, vmax=1,
                           extent=[0, grid_size[1], 0, grid_size[0]])
        plt.colorbar(im0, ax=ax[0], label='Allele Frequency')
        ax[0].set_title(f"Heat-Tolerance Allele\nYear {t * dt:.0f}")

        # Population density
        im1 = ax[1].imshow(N.T, cmap='viridis',
                           norm=plt.Normalize(0, K0),
                           extent=[0, grid_size[1], 0, grid_size[0]])
        plt.colorbar(im1, ax=ax[1], label='Population Density')
        ax[1].set_title("Population Density")

        # Environment
        im2 = ax[2].imshow(C.T, cmap='Greens',
                           extent=[0, grid_size[1], 0, grid_size[0]])
        plt.colorbar(im2, ax=ax[2], label='Habitat Conductance')
        ax[2].set_title("Habitat Quality")

        plt.tight_layout()
        plt.savefig(f'temperate_forest_simulation/frame_{frame:03d}.png', dpi=120)
        plt.close()
        frame += 1

print(f"Success! {frame} frames saved with:")
print("- Fixed urban corridor assignment")
print("- Stable population dynamics")
print("- Clear visualization labels")