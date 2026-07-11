"""
Quick plotting helper for the plain-text output written by the
Fortran solvers (results/*.dat).

Usage:
    python scripts/plot_results.py poisson2d results/poisson2d_custom_source.dat
    python scripts/plot_results.py transport1d results/transport1d_t*.dat

Requires: numpy, matplotlib  (pip install numpy matplotlib)
"""
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt


def load_2d(path):
    data = np.loadtxt(path)
    x, y, u = data[:, 0], data[:, 1], data[:, 2]
    nx = len(np.unique(x))
    ny = len(np.unique(y))
    X = x.reshape(ny, nx)
    Y = y.reshape(ny, nx)
    U = u.reshape(ny, nx)
    return X, Y, U


def plot_poisson2d(path):
    X, Y, U = load_2d(path)
    fig, ax = plt.subplots(figsize=(6, 5))
    cs = ax.contourf(X, Y, U, levels=40, cmap="RdBu_r")
    fig.colorbar(cs, ax=ax, label="u(x,y)")
    ax.set_title(path)
    ax.set_xlabel("x"); ax.set_ylabel("y")
    ax.set_aspect("equal")
    plt.tight_layout()
    plt.savefig(path.replace(".dat", ".png"), dpi=150)
    print("saved", path.replace(".dat", ".png"))


def plot_transport1d(paths):
    paths = sorted(paths)
    fig, ax = plt.subplots(figsize=(6, 4))
    for i, p in enumerate(paths):
        data = np.loadtxt(p)
        alpha = 0.15 + 0.85 * i / max(1, len(paths) - 1)
        ax.plot(data[:, 0], data[:, 1], color="C0", alpha=alpha)
    ax.set_xlabel("x"); ax.set_ylabel("u(x,t)")
    ax.set_title("1D transport: early (light) -> late (dark)")
    plt.tight_layout()
    plt.savefig("results/transport1d_overlay.png", dpi=150)
    print("saved results/transport1d_overlay.png")


def plot_transport2d_last(paths):
    path = sorted(paths)[-1]
    plot_poisson2d(path)  # same contour-plot machinery works fine


if __name__ == "__main__":
    kind = sys.argv[1]
    file_args = sys.argv[2:]
    expanded = []
    for f in file_args:
        expanded.extend(glob.glob(f))

    if kind == "poisson2d":
        plot_poisson2d(expanded[0])
    elif kind == "transport1d":
        plot_transport1d(expanded)
    elif kind == "transport2d":
        plot_transport2d_last(expanded)
    else:
        print("unknown kind:", kind)
