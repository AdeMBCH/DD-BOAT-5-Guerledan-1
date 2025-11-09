# Code généré par CHATGPT A partir du script de calibration pour le plot.

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from mpl_toolkits.mplot3d import Axes3D

DATA_PATH = "datacalibration13.txt"
OUT_BEFORE = "nuage_avant_xyz.png"
OUT_AFTER  = "nuage_apres_xyz.png"

def fit_ellipsoid_ls(X: np.ndarray):
    x, y, z = X[:, 0], X[:, 1], X[:, 2]
    D = np.column_stack([x*x, y*y, z*z, x*y, x*z, y*z, x, y, z, np.ones_like(x)])
    _, _, Vt = np.linalg.svd(D, full_matrices=False)
    p = Vt[-1, :]

    Ae = np.array([
        [p[0],     p[3]/2.0, p[4]/2.0],
        [p[3]/2.0, p[1],     p[5]/2.0],
        [p[4]/2.0, p[5]/2.0, p[2]]
    ], dtype=float)
    b = np.array([p[6], p[7], p[8]], dtype=float)
    c = float(p[9])

    def solve_once(sign=1.0):
        Ae_s = Ae * sign
        try:
            o = -0.5 * np.linalg.solve(Ae_s, b * sign)
        except np.linalg.LinAlgError:
            return None
        denom = float(o.T @ Ae_s @ o - c * sign)
        if abs(denom) < 1e-18:
            return None
        k = 1.0 / denom
        Ae_k = Ae_s * k
        try:
            L = np.linalg.cholesky(Ae_k)
        except np.linalg.LinAlgError:
            return None
        A = L.T
        return o, A, Ae_k

    out = solve_once(1.0) or solve_once(-1.0)
    if out is None:
        raise RuntimeError("Échec ajustement ellipsoïde.")
    return out

def calibrate_xyz(raw_xyz: np.ndarray):
    o, A, _ = fit_ellipsoid_ls(raw_xyz)
    Y = (raw_xyz - o) @ A.T
    return o, A, Y

def load_xyz(path: str) -> np.ndarray:
    df = pd.read_csv(Path(path), sep=None, engine="python")
    cols = ["xmag", "ymag", "zmag"]
    if all(c in df.columns for c in cols):
        M = df[cols].to_numpy(dtype=float)
    else:
        num = df.select_dtypes(include=[np.number])
        if num.shape[1] < 3:
            raise ValueError("3 colonnes numériques requises.")
        M = num.iloc[:, :3].to_numpy(dtype=float)
    return M

def set_axes_equal(ax):
    # Aspect 1:1:1 en 3D
    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()
    x_range = x_limits[1] - x_limits[0]
    y_range = y_limits[1] - y_limits[0]
    z_range = z_limits[1] - z_limits[0]
    max_range = max([x_range, y_range, z_range])
    x_middle = np.mean(x_limits)
    y_middle = np.mean(y_limits)
    z_middle = np.mean(z_limits)
    ax.set_xlim3d([x_middle - max_range/2, x_middle + max_range/2])
    ax.set_ylim3d([y_middle - max_range/2, y_middle + max_range/2])
    ax.set_zlim3d([z_middle - max_range/2, z_middle + max_range/2])

def main():
    M = load_xyz(DATA_PATH)
    o, A, Y = calibrate_xyz(M)

    # Figure 1: avant
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111, projection='3d')
    ax1.scatter(M[:,0], M[:,1], M[:,2], s=4)
    ax1.set_xlabel("X mag (raw)")
    ax1.set_ylabel("Y mag (raw)")
    ax1.set_zlabel("Z mag (raw)")
    ax1.set_title("Avant calibration (XYZ)")
    set_axes_equal(ax1)
    fig1.savefig(OUT_BEFORE, dpi=150, bbox_inches="tight")

    # Figure 2: après
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111, projection='3d')
    ax2.scatter(Y[:,0], Y[:,1], Y[:,2], s=4)
    # Sphère unité (wireframe)
    u = np.linspace(0, 2*np.pi, 60)
    v = np.linspace(0, np.pi, 30)
    xs = np.outer(np.cos(u), np.sin(v))
    ys = np.outer(np.sin(u), np.sin(v))
    zs = np.outer(np.ones_like(u), np.cos(v))
    ax2.plot_wireframe(xs, ys, zs, rstride=4, cstride=4, linewidth=0.4, alpha=0.4)

    ax2.set_xlabel("X mag calibré")
    ax2.set_ylabel("Y mag calibré")
    ax2.set_zlabel("Z mag calibré")
    ax2.set_title("Calibré + sphère unité")
    set_axes_equal(ax2)
    fig2.savefig(OUT_AFTER, dpi=150, bbox_inches="tight")

    nb, na = np.linalg.norm(M, axis=1), np.linalg.norm(Y, axis=1)
    print("Biais o:", o)
    print("Matrice A:\n", A)
    print("Norme moyenne avant:", nb.mean(), "écart-type:", nb.std())
    print("Norme moyenne après :", na.mean(), "écart-type:", na.std())
    print("Images:", OUT_BEFORE, "et", OUT_AFTER)

    plt.show()

if __name__ == "__main__":
    main()
