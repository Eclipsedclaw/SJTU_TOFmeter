#!/usr/bin/env python3
"""
pedestal_analysis.py

Reusable pedestal statistics + histogram + Gaussian-fit tool for CAEN
digitizer WaveDump output (wave0.txt ... waveN.txt, one sample per line).

Usage:
    python3 pedestal_analysis.py /path/to/folder
    python3 pedestal_analysis.py /path/to/folder --channels 16 --label "Slot 3"
    python3 pedestal_analysis.py /path/to/folder --out /path/to/outdir

Outputs (written into --out, default: same as input folder):
    pedestal_histograms.png   - 4x4 (or NxN) grid of histograms + Gaussian fits
    pedestal_summary.csv      - per-channel raw stats + fit parameters
"""

import argparse
import math
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def gaussian(x, amplitude, mean, sigma):
    return amplitude * np.exp(-0.5 * ((x - mean) / sigma) ** 2)


def fit_gaussian(vals, bins):
    """Fit a Gaussian to the histogram of vals. Returns dict with fit
    results, or None if the fit could not be performed (e.g. flat channel
    or fit failure)."""
    counts, edges = np.histogram(vals, bins=bins)
    centers = 0.5 * (edges[:-1] + edges[1:])

    if np.count_nonzero(counts) < 3:
        return None  # not enough distinct bins to fit meaningfully

    p0 = [counts.max(), float(np.mean(vals)), float(np.std(vals)) or 1.0]
    try:
        popt, pcov = curve_fit(gaussian, centers, counts, p0=p0, maxfev=10000)
        perr = np.sqrt(np.diag(pcov)) if pcov is not None else [np.nan] * 3
        amp, mean, sigma = popt
        sigma = abs(sigma)
        return {
            "amplitude": amp, "mean": mean, "sigma": sigma,
            "amplitude_err": perr[0], "mean_err": perr[1], "sigma_err": perr[2],
            "centers": centers, "counts": counts,
        }
    except Exception:
        return None


def analyze_channel(path):
    vals = np.loadtxt(path)
    if vals.ndim == 0:
        vals = np.array([float(vals)])
    return vals


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input_folder", help="Folder containing wave0.txt ... waveN.txt")
    ap.add_argument("--channels", type=int, default=16,
                     help="Number of channels to look for (default: 16)")
    ap.add_argument("--prefix", default="wave",
                     help="Filename prefix before the channel number (default: 'wave')")
    ap.add_argument("--out", default=None,
                     help="Output directory (default: same as input_folder)")
    ap.add_argument("--label", default=None,
                     help="Label for the plot title, e.g. 'Slot 3 (S/N 426)'")
    args = ap.parse_args()

    in_dir = args.input_folder
    out_dir = args.out or in_dir
    os.makedirs(out_dir, exist_ok=True)

    n = args.channels
    data = {}
    missing = []
    for ch in range(n):
        path = os.path.join(in_dir, f"{args.prefix}{ch}.txt")
        if not os.path.isfile(path):
            missing.append(ch)
            continue
        data[ch] = analyze_channel(path)

    if missing:
        print(f"WARNING: missing channel files, skipped: {missing}", file=sys.stderr)

    present_channels = sorted(data.keys())
    if not present_channels:
        print("ERROR: no channel files found.", file=sys.stderr)
        sys.exit(1)

    # ---- Compute stats + fits, print + collect for CSV ----
    header = (f"{'Ch':>3} {'N':>5} {'RawMean':>10} {'RawStd':>8} "
              f"{'FitMean':>10} {'FitSigma':>9} {'FitAmp':>9} {'Status':>10}")
    print(header)
    print("-" * len(header))

    rows = []
    for ch in present_channels:
        vals = data[ch]
        raw_mean = float(np.mean(vals))
        raw_std = float(np.std(vals))
        n_samp = len(vals)

        if raw_std == 0:
            status = "FLAT/DEAD"
            fit = None
        else:
            lo = int(math.floor(vals.min())) - 1
            hi = int(math.ceil(vals.max())) + 2
            bins = np.arange(lo, hi, 1)
            fit = fit_gaussian(vals, bins)
            status = "OK" if fit else "FIT_FAILED"

        if fit:
            print(f"{ch:>3} {n_samp:>5} {raw_mean:>10.2f} {raw_std:>8.3f} "
                  f"{fit['mean']:>10.2f} {fit['sigma']:>9.3f} {fit['amplitude']:>9.1f} {status:>10}")
        else:
            print(f"{ch:>3} {n_samp:>5} {raw_mean:>10.2f} {raw_std:>8.3f} "
                  f"{'--':>10} {'--':>9} {'--':>9} {status:>10}")

        rows.append({
            "channel": ch, "n_samples": n_samp,
            "raw_mean": raw_mean, "raw_std": raw_std,
            "fit_mean": fit["mean"] if fit else "",
            "fit_mean_err": fit["mean_err"] if fit else "",
            "fit_sigma": fit["sigma"] if fit else "",
            "fit_sigma_err": fit["sigma_err"] if fit else "",
            "fit_amplitude": fit["amplitude"] if fit else "",
            "status": status,
        })

    # ---- Write CSV ----
    csv_path = os.path.join(out_dir, "pedestal_summary.csv")
    with open(csv_path, "w") as f:
        cols = list(rows[0].keys())
        f.write(",".join(cols) + "\n")
        for r in rows:
            f.write(",".join(str(r[c]) for c in cols) + "\n")
    print(f"\nSaved summary CSV to {csv_path}")

    # ---- Plot grid ----
    ncols = 4
    nrows = math.ceil(n / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3 * nrows))
    axes = np.atleast_1d(axes).flatten()

    for ch in range(n):
        ax = axes[ch]
        if ch not in data:
            ax.set_title(f"Ch{ch}: NO DATA", fontsize=10, color="gray")
            ax.axis("off")
            continue

        vals = data[ch]
        raw_mean = float(np.mean(vals))
        raw_std = float(np.std(vals))

        if raw_std == 0:
            ax.hist(vals, bins=1, color="gray", alpha=0.7)
            ax.set_title(f"Ch{ch}: FLAT (mean={raw_mean:.0f})", fontsize=9, color="red")
        else:
            lo = int(math.floor(vals.min())) - 1
            hi = int(math.ceil(vals.max())) + 2
            bins = np.arange(lo, hi, 1)
            ax.hist(vals, bins=bins, color="#3b6fa0", alpha=0.8,
                    edgecolor="white", linewidth=0.3, label="data")

            fit = fit_gaussian(vals, bins)
            if fit:
                xx = np.linspace(lo, hi, 300)
                ax.plot(xx, gaussian(xx, fit["amplitude"], fit["mean"], fit["sigma"]),
                        color="red", linewidth=1.5, label="Gaussian fit")
                title = f"Ch{ch}: μ={fit['mean']:.1f}, σ={fit['sigma']:.2f}"
            else:
                ax.axvline(raw_mean, color="red", linestyle="--", linewidth=1)
                title = f"Ch{ch}: μ={raw_mean:.1f}, σ={raw_std:.2f} (no fit)"

            ax.set_title(title, fontsize=9)

        ax.set_xlabel("ADC counts", fontsize=7)
        ax.set_ylabel("Events", fontsize=7)
        ax.tick_params(labelsize=6)

    # hide any unused axes beyond n
    for extra in range(n, len(axes)):
        axes[extra].axis("off")

    suptitle = "Pedestal Histograms"
    if args.label:
        suptitle += f" — {args.label}"
    plt.suptitle(suptitle, fontsize=13, y=1.00)
    plt.tight_layout()

    png_path = os.path.join(out_dir, "pedestal_histograms.png")
    plt.savefig(png_path, dpi=150, bbox_inches="tight")
    print(f"Saved plot to {png_path}")


if __name__ == "__main__":
    main()