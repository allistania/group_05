# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import pandas as pd
import sys
import os

def plot_quality_histograms(csv_raw_file, hist_dat_file, output_prefix=None):

    df = pd.read_csv(csv_raw_file)
    ar = df['AR']
    skew = df['Skewness']

    fig, axes = plt.subplots(1, 3, figsize=(15, 4))

    axes[0].hist(ar, bins=30, edgecolor='black', alpha=0.7)
    axes[0].set_xlabel('Aspect Ratio')
    axes[0].set_ylabel('Count')
    axes[0].set_title('Aspect Ratio Distribution')
    axes[0].axvline(ar.mean(), color='red', linestyle='dashed', linewidth=1, label=f'Mean = {ar.mean():.2f}')
    axes[0].legend()

    axes[1].hist(skew, bins=30, edgecolor='black', alpha=0.7, color='orange')
    axes[1].set_xlabel('Skewness')
    axes[1].set_ylabel('Count')
    axes[1].set_title('Skewness Distribution')
    axes[1].axvline(skew.mean(), color='red', linestyle='dashed', linewidth=1, label=f'Mean = {skew.mean():.2f}')
    axes[1].legend()

    # Scatter plot
    axes[2].scatter(ar, skew, s=1, alpha=0.5)
    axes[2].set_xlabel('Aspect Ratio')
    axes[2].set_ylabel('Skewness')
    axes[2].set_title('Skewness vs Aspect Ratio')
    axes[2].grid(True, linestyle='--', alpha=0.5)

    plt.tight_layout()

    if output_prefix:
        fig.savefig(f"{output_prefix}_quality_plots.png", dpi=150)
        print(f"Plots saved to {output_prefix}_quality_plots.png")
    else:
        plt.show()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python plot_quality.py <path_to_quality_raw.csv> [output_prefix]")
        sys.exit(1)

    csv_file = sys.argv[1]
    prefix = sys.argv[2] if len(sys.argv) > 2 else os.path.splitext(csv_file)[0]

    hist_file = prefix + "_quality_hist.dat"  # можно не использовать

    plot_quality_histograms(csv_file, hist_file, prefix)