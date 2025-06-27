import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import os

# Set working directory to script location
os.chdir(os.path.dirname(os.path.abspath(__file__)))

# Use scientific plot style
plt.style.use('scientific')

# Load header by skipping lines starting with '#'
with open("./outputFULL.txt") as f:
    for line in f:
        if not line.startswith("#"):
            header = line.strip().split(",")
            break

# Load data
data = np.loadtxt("./outputFULL.txt", delimiter=",", skiprows=15)

# Get column indices
mD_idx = header.index("mD")
mS_idx = header.index("mS")
Omega_idx = header.index("Omegah2")
yChi_idx = header.index("YChi")
ThetaY_idx = header.index("ThetaY")

# Loop over yChi values
for yChi_val in np.arange(1.0, 2.1, 1.0):
    for thetaY_val in [0., 0.1, 0.2, 0.3, 0.4, 0.5]:  # Fixed ThetaY value as per your data
        # Also fix ThetaY value (e.g., 0.0 as in your data)
        filtered_data = data[
            np.isclose(data[:, yChi_idx], yChi_val) &
            np.isclose(data[:, ThetaY_idx], thetaY_val, atol=1e-2)
        ]

        if filtered_data.size == 0:
            print(f"No data for YChi = {yChi_val:.1f}, ThetaY = {thetaY_val:.2f}")
            continue

        # Extract relevant columns
        x = filtered_data[:, mD_idx]
        y = filtered_data[:, mS_idx]
        z = filtered_data[:, Omega_idx]

        # Create grid
        xi = np.unique(x)
        yi = np.unique(y)
        try:
            Xi, Yi = np.meshgrid(xi, yi)
            Zi = z.reshape(len(yi), len(xi))
        except ValueError:
            print(f"Skipping yChi = {yChi_val:.1f}, ThetaY = {thetaY_val:.2f}: irregular grid")
            continue

        # Set color normalization and levels
        norm = mcolors.LogNorm(vmin=1e-3, vmax=1.0, clip=False)
        levels = np.logspace(-3, 0, 11)  # 10 log-spaced contours

        # Plot
        plt.figure(figsize=(5, 4))
        cb = plt.contourf(Xi, Yi, Zi, levels=levels, cmap='viridis', norm=norm, extend='max')
        contour = plt.contour(Xi,Yi,Zi, levels=[0.12], colors='red', linewidths=1, linestyles='dashed')

        # Add colorbar with fixed ticks
        cbar = plt.colorbar(cb, label=r"$\Omega_\text{DM}h^2$")
        cbar.set_ticks([1e-3, 1e-2, 1e-1, 1e0])
        cbar.set_ticklabels([r"$10^{-3}$", r"$10^{-2}$", r"$10^{-1}$", r"$10^{0}$"])

        # Axis labels and title
        plt.xlabel(r"$m_D$ [GeV]")
        plt.ylabel(r"$m_S$ [GeV]")
        plt.xlim(1000,3000)
        plt.ylim(1000,3000)
        plt.title(rf"Relic density for $Y_{{\chi}} = {yChi_val:.1f}$ and $\tan{{\theta}} = {np.tan(thetaY_val):.2f}$")
        # Label the contour line
        if len(contour.allsegs[0]) > 0:
            plt.clabel(contour, fmt={0.12: r'$\Omega_\text{DM}h^2 = 0.12$'}, colors='red', fontsize=12, inline=True)
        else:
            print("No contour at level 0.12 to label.")
        plt.tight_layout()
        plt.savefig(f"figures/RelicDensity_Y={yChi_val:.1f}TanTheta={np.tan(thetaY_val):.1f}.pdf")
        plt.close()
