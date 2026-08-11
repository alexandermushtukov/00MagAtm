#!/usr/bin/env python3
import os
from io import BytesIO

import numpy as np
import matplotlib.pyplot as plt
import imageio.v2 as imageio


INPUT_FILE = "./res/res_B1e14_theta0_T1_40_30_12_iter_"
OUTPUT_GIF = "./res/res_B1e14_theta0_T1_40_30_12_iter_.gif"

# Global style
plt.rcParams.update({
    "font.size": 16,
    "axes.labelsize": 18,
    "axes.titlesize": 18,
    "xtick.labelsize": 15,
    "ytick.labelsize": 15,
    "legend.fontsize": 14,
})

# Scientific-looking colors
DARK_GRAY = "#3a3a3a"
ORANGE = "#e68613"
LIGHT_BLUE = "#5fa8d3"


def parse_tables(path):
    """
    Read the input file and extract numeric tables.

    Expected format:
    - comment/separator lines start with '#'
    - each data row has exactly 9 numeric columns
    - tables are separated by comment lines and/or blank lines

    Returns:
        list of numpy arrays with shape (N, 9)
    """
    tables = []
    current = []

    with open(path, "r", encoding="utf-8") as f:
        for raw_line in f:
            line = raw_line.strip()

            if not line:
                if current:
                    tables.append(np.array(current, dtype=float))
                    current = []
                continue

            if line.startswith("#"):
                if current:
                    tables.append(np.array(current, dtype=float))
                    current = []
                continue

            parts = line.split()
            if len(parts) != 9:
                if current:
                    tables.append(np.array(current, dtype=float))
                    current = []
                continue

            try:
                row = [float(x) for x in parts]
            except ValueError:
                if current:
                    tables.append(np.array(current, dtype=float))
                    current = []
                continue

            current.append(row)

    if current:
        tables.append(np.array(current, dtype=float))

    return tables


def padded_limits(values, pad_fraction=0.05):
    """
    Build slightly padded linear-axis limits.
    """
    vmin = np.min(values)
    vmax = np.max(values)

    if np.isclose(vmin, vmax):
        delta = 1.0 if np.isclose(vmin, 0.0) else 0.05 * abs(vmin)
        return vmin - delta, vmax + delta

    pad = pad_fraction * (vmax - vmin)
    return vmin - pad, vmax + pad


def compute_global_limits(tables):
    """
    Compute fixed y-limits for each panel over all tables,
    so the vertical axes do not change from frame to frame.
    """
    panel1_all = []
    panel2_all = []
    panel3_all = []

    for table in tables:
        y1 = table[:, 3]  # column 4

        with np.errstate(divide="ignore", invalid="ignore"):
            y2 = table[:, 5] / table[:, 6]  # column 6 / column 7

        y3a = table[:, 5]  # column 6
        y3b = table[:, 7]  # column 8
        y3c = table[:, 8]  # column 9

        panel1_all.append(y1[np.isfinite(y1)])
        panel2_all.append(y2[np.isfinite(y2)])
        panel3_all.append(y3a[np.isfinite(y3a)])
        panel3_all.append(y3b[np.isfinite(y3b)])
        panel3_all.append(y3c[np.isfinite(y3c)])

    panel1_all = np.concatenate(panel1_all)
    panel2_all = np.concatenate(panel2_all)
    panel3_all = np.concatenate(panel3_all)

    return {
        "panel1": padded_limits(panel1_all),
        "panel2": padded_limits(panel2_all),
        "panel3": padded_limits(panel3_all),
    }


def add_slide_number(ax, slide_index):
    """
    Draw slide number in the top-left corner of the top panel only.
    """
    ax.text(
        0.02, 0.95,
        f"{slide_index}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=18,
        fontweight="bold",
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.7, pad=2.0),
    )


def make_frame(table, frame_index, total_frames, limits):
    """
    Build one GIF frame with three vertically stacked panels.
    """
    # Human-readable column numbering:
    # 1:i, 2:m, 3:rho, 4:T_keV, 5:dT, 6:F, 7:F_target, 8:J1, 9:J2
    x = table[:, 1]          # column 2
    y1 = table[:, 3]         # column 4

    with np.errstate(divide="ignore", invalid="ignore"):
        y2 = table[:, 5] / table[:, 6]   # column 6 / column 7

    y3a = table[:, 5]        # column 6
    y3b = table[:, 7]        # column 8
    y3c = table[:, 8]        # column 9

    fig, axes = plt.subplots(
        3, 1,
        figsize=(12, 11),
        sharex=True,
        constrained_layout=True
    )

    ax1, ax2, ax3 = axes

    # Top panel
    ax1.plot(x, y1, lw=2.5, color=DARK_GRAY)
    ax1.set_xscale("log")
    ax1.set_ylabel("T_keV")
    ax1.set_ylim(*limits["panel1"])
    ax1.grid(True, which="both", alpha=0.3)
    add_slide_number(ax1, frame_index + 1)

    # Middle panel
    ax2.plot(x, y2, lw=2.5, color=DARK_GRAY)
    ax2.axhline(1.0, color="black", linestyle=":", linewidth=1.8, alpha=0.8)
    ax2.set_xscale("log")
    ax2.set_ylabel(r"$F_i/F_{target}$")
    ax2.set_ylim(*limits["panel2"])
    ax2.grid(True, which="both", alpha=0.3)

    # Bottom panel
    ax3.plot(x, y3a, lw=2.5, color=DARK_GRAY, label="total flux")
    ax3.plot(x, y3b, lw=2.5, color=ORANGE, label="pol 1")
    ax3.plot(x, y3c, lw=2.5, color=LIGHT_BLUE, label="pol 2")
    ax3.set_xscale("log")
    ax3.set_ylabel("F")
    ax3.set_xlabel(r"m [g/cm^2]")
    ax3.set_ylim(*limits["panel3"])
    ax3.grid(True, which="both", alpha=0.3)
    ax3.legend(loc="best")

    # Save figure to memory
    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=120)
    plt.close(fig)
    buf.seek(0)

    return imageio.imread(buf)


def main():
    os.makedirs("./res", exist_ok=True)

    tables = parse_tables(INPUT_FILE)
    if not tables:
        raise RuntimeError(f"No tables were found in {INPUT_FILE}")

    limits = compute_global_limits(tables)

    frames = []
    total = len(tables)

    for i, table in enumerate(tables):
        if table.shape[1] != 9:
            print(f"Skipping table {i + 1}: invalid shape {table.shape}")
            continue

        frame = make_frame(table, i, total, limits)
        frames.append(frame)
        print(f"Prepared frame {i + 1}/{total}")

    if not frames:
        raise RuntimeError("No valid frames were generated")

    imageio.mimsave(OUTPUT_GIF, frames, duration=0.7, loop=0)
    print(f"GIF saved to: {OUTPUT_GIF}")


if __name__ == "__main__":
    main()
