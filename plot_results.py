#!/usr/bin/env python3

import csv
import os
from pathlib import Path


ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT / "output"
PLOTS = OUTPUT / "plots"
MPLCONFIG = OUTPUT / ".mplconfig"
CACHE_DIR = OUTPUT / ".cache"

MPLCONFIG.mkdir(parents=True, exist_ok=True)
CACHE_DIR.mkdir(parents=True, exist_ok=True)
os.environ["MPLCONFIGDIR"] = str(MPLCONFIG)
os.environ["XDG_CACHE_HOME"] = str(CACHE_DIR)

import matplotlib.pyplot as plt


def read_scalar_history(filename):
    rows = []
    with (OUTPUT / filename).open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            rows.append(
                {
                    "iteration": int(row["iteration"]),
                    "estimate": float(row["estimate"]),
                    "residual": float(row["residual"]),
                    "error": float(row["error"]),
                    "rate": None if row["rate"] == "" else float(row["rate"]),
                    "phase": row["phase"],
                }
            )
    return rows


def plot_scalar_method(filename, title, output_name, switch_iteration=None):
    rows = read_scalar_history(filename)
    iterations = [row["iteration"] for row in rows]
    estimates = [row["estimate"] for row in rows]
    errors = [abs(row["error"]) for row in rows]
    rate_iterations = [row["iteration"] for row in rows if row["rate"] is not None]
    rates = [row["rate"] for row in rows if row["rate"] is not None]

    figure, axes = plt.subplots(3, 1, figsize=(8, 10), constrained_layout=True)
    figure.suptitle(title, fontsize=14, fontweight="bold")

    axes[0].plot(iterations, estimates, marker="o", linewidth=1.5, markersize=4)
    axes[0].set_ylabel("h (W/m^2-K)")
    axes[0].set_xlabel("Iteration")
    axes[0].grid(True, alpha=0.3)

    axes[1].plot(iterations, errors, marker="o", linewidth=1.5, markersize=4, color="darkorange")
    axes[1].set_ylabel("Absolute error")
    axes[1].set_xlabel("Iteration")
    axes[1].set_yscale("symlog", linthresh=1e-12)
    axes[1].grid(True, alpha=0.3)

    axes[2].plot(rate_iterations, rates, marker="o", linewidth=1.5, markersize=4, color="forestgreen")
    axes[2].set_ylabel("Rate")
    axes[2].set_xlabel("Iteration")
    axes[2].set_yscale("symlog", linthresh=1e-12)
    axes[2].grid(True, alpha=0.3)

    if switch_iteration is not None:
        for axis in axes:
            axis.axvline(switch_iteration, color="crimson", linestyle="--", linewidth=1.0)
            axis.text(
                switch_iteration + 0.2,
                axis.get_ylim()[1] * 0.85,
                "switch",
                color="crimson",
                fontsize=9,
            )

    figure.savefig(PLOTS / output_name, dpi=200)
    plt.close(figure)


def main():
    PLOTS.mkdir(parents=True, exist_ok=True)
    plot_scalar_method("part1_bisection.csv", "Part 1(a): Bisection Method", "bisection_history.png")
    plot_scalar_method("part1_newton.csv", "Part 1(b): Newton Method", "newton_history.png")
    plot_scalar_method(
        "part1_hybrid.csv",
        "Part 1(c): Hybrid Bisection + Newton Method",
        "hybrid_history.png",
        switch_iteration=10,
    )


if __name__ == "__main__":
    main()
