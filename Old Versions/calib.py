#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

DATA_FILE = "datacalibration4.txt"
G0 = 9.80

def load_mag(path=DATA_FILE):
    rows = []
    with open(path, "r", encoding="ascii") as f:
        _ = f.readline()  # en-tête
        for line in f:
            s = line.strip()
            if not s:
                continue
            s = s.replace("(", "").replace(")", "")
            p = s.split(",")
            if len(p) < 9:
                continue
            mx, my, mz = float(p[0]), float(p[1]), float(p[2])  # colonnes mag
            rows.append((mx, my, mz))
    X = np.array(rows, float)
    if X.size == 0:
        raise RuntimeError("Aucune donnée magnétomètre.")
    return X

def M(X):
    x, y, z = X[:,0], X[:,1], X[:,2]
    return np.column_stack([x*x, y*y, z*z, x*y, x*z, y*z, x, y, z, np.ones_like(x)])

def RMS(M):
    # plus robuste que (M^T M)^{-1} M^T
    return np.linalg.lstsq(M, np.ones(M.shape[0]), rcond=None)[0]

def pinv(Q):
    # pseudo-inverse
    w, V = np.linalg.eigh(0.5*(Q+Q.T))
    tol = max(1e-12, 1e-6*np.max(np.abs(w)))
    w_clip = np.where(w > tol, w, tol)
    Qinv = V @ np.diag(1.0 / w_clip) @ V.T
    return Qinv

def compute_A_and_minus_b(X, scale=1.0):
    p = RMS(M(X))
    p1,p2,p3,p4,p5,p6,p7,p8,p9,p10 = p
    Q = np.array([[p1, 0.5*p4, 0.5*p5],
                  [0.5*p4, p2,  0.5*p6],
                  [0.5*p5, 0.5*p6, p3]], float)
    ell = np.array([p7, p8, p9], float)
    Q = 0.5*(Q+Q.T)

    Qinv = pinv(Q)
    b = -0.5 * (Qinv @ ell)
    
    kappa = 1.0 - p10 + 0.25 * (ell @ (Qinv @ ell))
    if kappa <= 0.0:
        Q  = -Q
        ell = -ell
        p10 = -p10
        Qinv = pinv(Q)
        b = -0.5 * (Qinv @ ell)
        kappa = 1.0 - p10 + 0.25 * (ell @ (Qinv @ ell))

    Q_unit = Q / kappa

    # A = scale * sqrt(Q_unit)
    w, V = np.linalg.eigh(Q_unit)
    w = np.maximum(w, 1e-12)
    A = scale * (V @ np.diag(np.sqrt(w)) @ V.T)

    return A, -b

def main():
    mag = load_mag()
    A_mag, bneg_mag = compute_A_and_minus_b(mag, scale=1.0)

    print("A_mag:\n", A_mag)
    print("bneg_mag:", bneg_mag)


if __name__ == "__main__":
    main()
