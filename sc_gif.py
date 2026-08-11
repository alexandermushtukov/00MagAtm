#!/usr/bin/env python3

import os
from io import BytesIO

import numpy as np
import matplotlib.pyplot as plt
import imageio.v2 as imageio


INPUT_FILE = "./res/res_B1e14_theta05_T1_40_30_12_iter_"
OUTPUT_GIF = "./res/res_B1e14_theta05_T1_40_30_12_iter_.gif"


# ============================================================
# Global plot style
# ============================================================

plt.rcParams.update({
    "font.size": 16,
    "axes.labelsize": 18,
    "axes.titlesize": 18,
    "xtick.labelsize": 15,
    "ytick.labelsize": 15,
    "legend.fontsize": 14,
})


# ============================================================
# Colors
# ============================================================

DARK_GRAY = "#3a3a3a"
MEDIUM_GRAY = "#777777"
LIGHT_GRAY = "#b0b0b0"

ORANGE = "#e68613"
LIGHT_BLUE = "#5fa8d3"

CURRENT_POINT = "#111111"


# ============================================================
# Read iteration file
# ============================================================

def parse_iterations(path):
    """
    Read atmospheric structures and convergence metrics.

    Expected structure for every iteration:

        # atmosphere structure, iter = ...
        # i, m, rho, T_keV, dT, F, F_target, flux_mode1, flux_mode2
        <numeric table with 9 columns>
        ...
        # : delta_surf, eps_Flog, eps_Fmax, eps_rough
        # <delta_surf> <eps_Flog> <eps_Fmax> <eps_rough>

    Returns
    -------
    iterations : list of dictionaries

        iterations[i]["table"]
            numpy array with shape (N, 9)

        iterations[i]["metrics"]
            numpy array:
            [delta_surf, eps_Flog, eps_Fmax, eps_rough]
    """

    iterations = []

    current_table = []
    current_iteration = None

    waiting_for_metrics = False

    with open(path, "r", encoding="utf-8") as f:

        for raw_line in f:

            line = raw_line.strip()

            if not line:
                continue

            # ------------------------------------------------
            # Start of a new iteration
            # ------------------------------------------------

            if line.startswith("# atmosphere structure"):

                current_table = []
                current_iteration = None

                try:
                    current_iteration = int(line.split("=")[1])
                except Exception:
                    current_iteration = len(iterations)

                continue

            # ------------------------------------------------
            # Header announcing convergence metrics
            # ------------------------------------------------

            if line.startswith("# :"):

                waiting_for_metrics = True
                continue

            # ------------------------------------------------
            # Read convergence metrics
            # ------------------------------------------------

            if waiting_for_metrics and line.startswith("#"):

                parts = line[1:].split()

                if len(parts) >= 4:

                    try:
                        metrics = np.array(
                            [float(parts[i]) for i in range(4)],
                            dtype=float
                        )

                        if current_table:

                            iterations.append({
                                "iteration": current_iteration,
                                "table": np.array(
                                    current_table,
                                    dtype=float
                                ),
                                "metrics": metrics
                            })

                    except ValueError:
                        pass

                waiting_for_metrics = False
                current_table = []

                continue

            # ------------------------------------------------
            # Other comment lines
            # ------------------------------------------------

            if line.startswith("#"):
                continue

            # ------------------------------------------------
            # Atmospheric table row
            # ------------------------------------------------

            parts = line.split()

            if len(parts) == 9:

                try:
                    row = [float(x) for x in parts]
                    current_table.append(row)

                except ValueError:
                    pass

    return iterations


# ============================================================
# Axis helpers
# ============================================================

def padded_limits(values, pad_fraction=0.05):
    """
    Return linear-axis limits with a small padding.
    """

    values = np.asarray(values)
    values = values[np.isfinite(values)]

    if len(values) == 0:
        return 0.0, 1.0

    vmin = np.min(values)
    vmax = np.max(values)

    if np.isclose(vmin, vmax):

        delta = (
            1.0
            if np.isclose(vmin, 0.0)
            else 0.05 * abs(vmin)
        )

        return vmin - delta, vmax + delta

    pad = pad_fraction * (vmax - vmin)

    return vmin - pad, vmax + pad


def padded_limits_include_zero(values, pad_fraction=0.08):
    """
    Same as padded_limits(), but make sure zero is visible.
    Useful for signed convergence quantities such as delta_surf.
    """

    values = np.asarray(values)
    values = values[np.isfinite(values)]

    if len(values) == 0:
        return -1.0, 1.0

    vmin = min(np.min(values), 0.0)
    vmax = max(np.max(values), 0.0)

    if np.isclose(vmin, vmax):
        return -1.0, 1.0

    pad = pad_fraction * (vmax - vmin)

    return vmin - pad, vmax + pad


# ============================================================
# Global axis limits
# ============================================================

def compute_global_limits(iterations):
    """
    Compute fixed limits for all panels.

    Fixed limits prevent axes from changing during GIF animation.
    """

    panel1_all = []
    panel2_all = []
    panel3_all = []

    metrics_all = []

    for item in iterations:

        table = item["table"]

        # Temperature
        y1 = table[:, 3]

        # Flux / target flux
        with np.errstate(divide="ignore", invalid="ignore"):
            y2 = table[:, 5] / table[:, 6]

        # Fluxes
        y3a = table[:, 5]
        y3b = table[:, 7]
        y3c = table[:, 8]

        panel1_all.append(y1[np.isfinite(y1)])
        panel2_all.append(y2[np.isfinite(y2)])

        panel3_all.append(y3a[np.isfinite(y3a)])
        panel3_all.append(y3b[np.isfinite(y3b)])
        panel3_all.append(y3c[np.isfinite(y3c)])

        metrics_all.append(item["metrics"])

    panel1_all = np.concatenate(panel1_all)
    panel2_all = np.concatenate(panel2_all)
    panel3_all = np.concatenate(panel3_all)

    metrics_all = np.array(metrics_all)

    delta_surf = metrics_all[:, 0]
    eps_Flog   = metrics_all[:, 1]
    eps_Fmax   = metrics_all[:, 2]
    eps_rough  = metrics_all[:, 3]

    return {

        "panel1": padded_limits(panel1_all),
        "panel2": padded_limits(panel2_all),
        "panel3": padded_limits(panel3_all),

        "delta_surf":
            padded_limits_include_zero(delta_surf),

        "flux_errors":
            padded_limits(
                np.concatenate([eps_Flog, eps_Fmax])
            ),

        "eps_rough":
            padded_limits(eps_rough),
    }


# ============================================================
# Iteration number
# ============================================================

def add_slide_number(ax, iteration_number):
    """
    Put the current iteration number in the top-left panel.
    """

    ax.text(
        0.02,
        0.95,
        f"iteration {iteration_number}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=17,
        fontweight="bold",
        bbox=dict(
            facecolor="white",
            edgecolor="none",
            alpha=0.75,
            pad=2.5
        ),
    )


# ============================================================
# Plot one quality-metric panel
# ============================================================

def plot_metric_panel(
    ax,
    iteration_numbers,
    values,
    current_index,
    ylabel,
    ylim,
    linestyle="-",
    show_zero=False
):
    """
    Draw one convergence metric versus iteration number.
    """

    ax.plot(
        iteration_numbers,
        values,
        color=MEDIUM_GRAY,
        linewidth=2.0,
        linestyle=linestyle
    )

    # Current position
    ax.plot(
        iteration_numbers[current_index],
        values[current_index],
        marker="o",
        markersize=9,
        color=CURRENT_POINT,
        zorder=10
    )

    if show_zero:
        ax.axhline(
            0.0,
            color=LIGHT_GRAY,
            linestyle=":",
            linewidth=1.5
        )

    ax.set_ylabel(ylabel)
    ax.set_ylim(*ylim)

    ax.grid(
        True,
        alpha=0.3
    )


# ============================================================
# Build one animation frame
# ============================================================

def make_frame(
    iterations,
    frame_index,
    limits
):
    """
    Build one GIF frame.

    Left column:
        1. Temperature profile
        2. F / F_target
        3. Total and polarization fluxes

    Right column:
        1. delta_surf
        2. eps_Flog and eps_Fmax
        3. eps_rough
    """

    item = iterations[frame_index]

    table = item["table"]

    iteration_number = item["iteration"]

    # --------------------------------------------------------
    # Atmospheric quantities
    # --------------------------------------------------------

    x = table[:, 1]

    # Temperature
    y1 = table[:, 3]

    # F / F_target
    with np.errstate(divide="ignore", invalid="ignore"):
        y2 = table[:, 5] / table[:, 6]

    # Flux components
    y3a = table[:, 5]
    y3b = table[:, 7]
    y3c = table[:, 8]

    # --------------------------------------------------------
    # Convergence history
    # --------------------------------------------------------

    iteration_numbers = np.array([
        x["iteration"] for x in iterations
    ])

    metrics = np.array([
        x["metrics"] for x in iterations
    ])

    delta_surf = metrics[:, 0]
    eps_Flog   = metrics[:, 1]
    eps_Fmax   = metrics[:, 2]
    eps_rough  = metrics[:, 3]

    # --------------------------------------------------------
    # Figure layout
    # --------------------------------------------------------

    fig = plt.figure(
        figsize=(17, 11),
        constrained_layout=True
    )

    gs = fig.add_gridspec(
        3,
        2,
        width_ratios=[3.1, 1.45]
    )

    # Left column
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0], sharex=ax1)
    ax3 = fig.add_subplot(gs[2, 0], sharex=ax1)

    # Right column
    ax4 = fig.add_subplot(gs[0, 1])
    ax5 = fig.add_subplot(gs[1, 1], sharex=ax4)
    ax6 = fig.add_subplot(gs[2, 1], sharex=ax4)

    # ========================================================
    # LEFT: atmospheric structure
    # ========================================================

    # --------------------------------------------------------
    # Panel 1: temperature
    # --------------------------------------------------------

    ax1.plot(
        x,
        y1,
        linewidth=2.5,
        color=DARK_GRAY
    )

    ax1.set_xscale("log")

    ax1.set_ylabel(r"$T_{\rm keV}$")

    ax1.set_ylim(
        *limits["panel1"]
    )

    ax1.grid(
        True,
        which="both",
        alpha=0.3
    )

    add_slide_number(
        ax1,
        iteration_number
    )

    # --------------------------------------------------------
    # Panel 2: normalized flux
    # --------------------------------------------------------

    ax2.plot(
        x,
        y2,
        linewidth=2.5,
        color=DARK_GRAY
    )

    ax2.axhline(
        1.0,
        color="black",
        linestyle=":",
        linewidth=1.8,
        alpha=0.8
    )

    ax2.set_xscale("log")

    ax2.set_ylabel(
        r"$F/F_{\rm target}$"
    )

    ax2.set_ylim(
        *limits["panel2"]
    )

    ax2.grid(
        True,
        which="both",
        alpha=0.3
    )

    # --------------------------------------------------------
    # Panel 3: flux components
    # --------------------------------------------------------

    ax3.plot(
        x,
        y3a,
        linewidth=2.5,
        color=DARK_GRAY,
        label="total flux"
    )

    ax3.plot(
        x,
        y3b,
        linewidth=2.5,
        color=ORANGE,
        label="pol 1"
    )

    ax3.plot(
        x,
        y3c,
        linewidth=2.5,
        color=LIGHT_BLUE,
        label="pol 2"
    )

    ax3.set_xscale("log")

    ax3.set_ylabel("F")

    ax3.set_xlabel(
        r"$m\ [{\rm g\,cm^{-2}}]$"
    )

    ax3.set_ylim(
        *limits["panel3"]
    )

    ax3.grid(
        True,
        which="both",
        alpha=0.3
    )

    # Fixed legend position
    ax3.legend(
        loc="lower left"
    )

    # Hide x tick labels in top two left panels
    plt.setp(
        ax1.get_xticklabels(),
        visible=False
    )

    plt.setp(
        ax2.get_xticklabels(),
        visible=False
    )

    # ========================================================
    # RIGHT: convergence / quality metrics
    # ========================================================

    # --------------------------------------------------------
    # Panel 4: delta_surf
    # --------------------------------------------------------

    plot_metric_panel(
        ax4,
        iteration_numbers,
        delta_surf,
        frame_index,
        r"$\delta_{\rm surf}$",
        limits["delta_surf"],
        show_zero=True
    )

    ax4.set_title(
        "Iteration quality"
    )

    # --------------------------------------------------------
    # Panel 5: eps_Flog + eps_Fmax
    # --------------------------------------------------------

    ax5.plot(
        iteration_numbers,
        eps_Flog,
        color=MEDIUM_GRAY,
        linewidth=2.0,
        linestyle="-",
        label=r"$\epsilon_{F,\log}$"
    )

    ax5.plot(
        iteration_numbers,
        eps_Fmax,
        color=MEDIUM_GRAY,
        linewidth=2.0,
        linestyle="--",
        label=r"$\epsilon_{F,\max}$"
    )

    # Current points
    ax5.plot(
        iteration_numbers[frame_index],
        eps_Flog[frame_index],
        marker="o",
        markersize=9,
        color=CURRENT_POINT,
        zorder=10
    )

    ax5.plot(
        iteration_numbers[frame_index],
        eps_Fmax[frame_index],
        marker="s",
        markersize=8,
        color=CURRENT_POINT,
        zorder=10
    )

    ax5.set_ylabel(
        r"$\epsilon_F$"
    )

    ax5.set_ylim(
        *limits["flux_errors"]
    )

    ax5.grid(
        True,
        alpha=0.3
    )

    ax5.legend(
        loc="upper right",
        fontsize=12
    )

    # --------------------------------------------------------
    # Panel 6: roughness
    # --------------------------------------------------------

    plot_metric_panel(
        ax6,
        iteration_numbers,
        eps_rough,
        frame_index,
        r"$\epsilon_{\rm rough}$",
        limits["eps_rough"]
    )

    ax6.set_xlabel(
        "Iteration"
    )

    # --------------------------------------------------------
    # Fixed x-range for convergence panels
    # --------------------------------------------------------

    if len(iteration_numbers) > 1:

        xmin = np.min(iteration_numbers)
        xmax = np.max(iteration_numbers)

        dx = xmax - xmin

        ax4.set_xlim(
            xmin - 0.03 * dx,
            xmax + 0.03 * dx
        )

    else:

        ax4.set_xlim(
            iteration_numbers[0] - 1,
            iteration_numbers[0] + 1
        )

    # Hide iteration labels on upper quality panels
    plt.setp(
        ax4.get_xticklabels(),
        visible=False
    )

    plt.setp(
        ax5.get_xticklabels(),
        visible=False
    )

    # --------------------------------------------------------
    # Save frame to memory
    # --------------------------------------------------------

    buf = BytesIO()

    fig.savefig(
        buf,
        format="png",
        dpi=120
    )

    plt.close(fig)

    buf.seek(0)

    return imageio.imread(buf)


# ============================================================
# Main
# ============================================================

def main():

    os.makedirs(
        "./res",
        exist_ok=True
    )

    # --------------------------------------------------------
    # Read data
    # --------------------------------------------------------

    iterations = parse_iterations(
        INPUT_FILE
    )

    if not iterations:

        raise RuntimeError(
            f"No iterations were found in {INPUT_FILE}"
        )

    print(
        f"Found {len(iterations)} iterations"
    )

    # --------------------------------------------------------
    # Compute global plot limits
    # --------------------------------------------------------

    limits = compute_global_limits(
        iterations
    )

    # --------------------------------------------------------
    # Generate frames
    # --------------------------------------------------------

    frames = []

    total = len(iterations)

    for i in range(total):

        table = iterations[i]["table"]

        if table.shape[1] != 9:

            print(
                f"Skipping iteration {i}: "
                f"invalid table shape {table.shape}"
            )

            continue

        frame = make_frame(
            iterations,
            i,
            limits
        )

        frames.append(
            frame
        )

        print(
            f"Prepared frame {i + 1}/{total}"
        )

    if not frames:

        raise RuntimeError(
            "No valid frames were generated"
        )

    # --------------------------------------------------------
    # Save GIF
    # --------------------------------------------------------

    imageio.mimsave(
        OUTPUT_GIF,
        frames,
        duration=0.7,
        loop=0
    )

    print(
        f"GIF saved to: {OUTPUT_GIF}"
    )


if __name__ == "__main__":
    main()
