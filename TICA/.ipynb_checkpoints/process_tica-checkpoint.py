#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.collections import LineCollection
import seaborn as sns
import matplotlib
import sys

# Set plot context for smaller font sizes
sns.set_context('talk', font_scale=1.1)  # 'talk' context is smaller than 'poster', and font_scale adjusts text size

# Check command-line arguments
if len(sys.argv) != 8:
    print('Usage: SCRIPT [H_path] [xedges_path] [yedges_path] [rmsd_path] [X_TICA1] [X_TICA2] [output_path]')
    sys.exit(2)

H_path = sys.argv[1]
xedges_path = sys.argv[2]
yedges_path = sys.argv[3]
rmsd_path = sys.argv[4]
X_TICA1_path = sys.argv[5]
X_TICA2_path = sys.argv[6]
output_path = sys.argv[7]

def plot_energy_landscape_with_rmsd(H_path, xedges_path, yedges_path, rmsd_path, X_TICA1_path, X_TICA2_path, output_path):
    # Load energy landscape data
    H_ML = np.genfromtxt(H_path, delimiter=' ')
    if H_ML.ndim == 1:
        H_ML = H_ML.reshape((50, 50))
    energy = -np.log(H_ML + 1e-16)
    
    # Prepare meshgrid for the energy landscape
    xedges_ML = np.loadtxt(xedges_path)
    yedges_ML = np.loadtxt(yedges_path)
    X, Y = np.meshgrid(xedges_ML[:-1], yedges_ML[:-1])
    
    # Load TICA vectors
    TICA1 = np.loadtxt(X_TICA1_path)
    TICA2 = np.loadtxt(X_TICA2_path)

    # Load RMSD values
    rmsd_values = np.loadtxt(rmsd_path)
    norm = Normalize(vmin=min(rmsd_values), vmax=max(rmsd_values))
    cmap = matplotlib.cm.get_cmap('rainbow')

    # Set up the figure with two subplots, with a larger figure size
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10), gridspec_kw={'height_ratios': [3, 1]})
    
    # Plot the energy landscape as a contour plot
    levels = np.linspace(0,
                     12,   # 0 → 最高能
                     13)            # 想要多少层就写多少
    contour = ax1.contourf(X, Y, energy.T, levels=levels, cmap='viridis')
    clb = fig.colorbar(contour, ax=ax1)
    clb.set_label('Free Energy (kcal/mol)')
    ax1.set_xlabel('TICA1')
    ax1.set_ylabel('TICA2')
    #ax1.set_title('Free Energy Landscape of ML')
    
    # Use the full range of xedges_ML and yedges_ML for axis limits to fit the energy landscape
    ax1.set_xlim([xedges_ML[0], xedges_ML[-1]])
    ax1.set_ylim([yedges_ML[0], yedges_ML[-1]])

    # Reduce the number of ticks on the y-axis
    ax1.set_yticks(np.linspace(yedges_ML[0], yedges_ML[-1], 5))
    ax1.set_xticks(np.linspace(xedges_ML[0], xedges_ML[-1], 5))# Adjust 5 to control the number of ticks
    ax1.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:.5f}'))
    ax1.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: f'{y:.5f}'))
    
    # Divide TICA data into segments and plot each with a gradient color
    skip = 1000  # Adjust to control segment length
    segments_x = np.stack((TICA1[::skip][:-1], TICA1[::skip][1:]))
    segments_y = np.stack((TICA2[::skip][:-1], TICA2[::skip][1:]))
    nsegs = segments_x.shape[1]
    
    for i in range(nsegs):
        ax1.plot(segments_x[:, i], segments_y[:, i], color=cmap(i / nsegs), linewidth=2)

    # Calculate time in microseconds, with each step as 0.01 μs
    time = np.arange(len(rmsd_values)) * 10 ** (-5)  # or simply `0.01` for readability

    # Plot RMSD over time with gradient coloring
    for i in range(nsegs):
        if i == nsegs - 1:
            ax2.plot(time[i * skip:], rmsd_values[i * skip:], color=cmap(i / nsegs))
        else:
            ax2.plot(time[i * skip:(i + 1) * skip + 1], rmsd_values[i * skip:(i + 1) * skip + 1], color=cmap(i / nsegs))
    
    ax2.set_xlabel('Time [μs]')  # Time unit in microseconds
    ax2.set_ylabel('RMSD [Å]')
    ax2.set_xlim([0, time[-1]])
    ax2.grid(visible=False)
    #ax2.set_title("RMSD Trajectory")

    # Print the ranges of TICA and energy landscape data
    print(f"TICA1 range: {TICA1.min()}, {TICA1.max()}")
    print(f"TICA2 range: {TICA2.min()}, {TICA2.max()}")
    print(f"xedges_ML range: {xedges_ML[0]}, {xedges_ML[-1]}")
    print(f"yedges_ML range: {yedges_ML[0]}, {yedges_ML[-1]}")
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight',format = 'pdf')
    print(f"Plot saved to {output_path}")

# Run the function with the provided paths
plot_energy_landscape_with_rmsd(H_path, xedges_path, yedges_path, rmsd_path, X_TICA1_path, X_TICA2_path, output_path)
sns.reset_orig()