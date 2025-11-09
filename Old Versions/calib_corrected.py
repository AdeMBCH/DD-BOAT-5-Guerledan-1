# calibrate_magnetometer.py
# Ellipsoid fitting LS -> bias o and gain matrix A (3D), then y = A @ (x - o)
# Méthode conforme aux équations linéaires d'ellipsoïde du papier.

import numpy as np
import pandas as pd
from pathlib import Path

def fit_ellipsoid_ls(X: np.ndarray):
    """
    X: Nx3 array of raw magnetometer samples [x, y, z]
    Returns:
      o  : 3x1 bias (hard-iron)
      A  : 3x3 calibration matrix such that y = A @ (x - o) and ||y|| = 1
      Ae : 3x3 symmetric matrix with A^T A = Ae
    """
    x, y, z = X[:, 0], X[:, 1], X[:, 2]

    # Design matrix D with d(x) = [x^2, y^2, z^2, x*y, x*z, y*z, x, y, z, 1]^T
    D = np.column_stack([x*x, y*y, z*z, x*y, x*z, y*z, x, y, z, np.ones_like(x)])

    # Solve D p = 0 via SVD (p = right singular vector for smallest singular value)
    U, S, Vt = np.linalg.svd(D, full_matrices=False)
    p = Vt[-1, :]  # 10 params

    def unpack(pvec):
        Ae = np.array([
            [pvec[0],       pvec[3]/2.0, pvec[4]/2.0],
            [pvec[3]/2.0,   pvec[1],     pvec[5]/2.0],
            [pvec[4]/2.0,   pvec[5]/2.0, pvec[2]]
        ], dtype=float)
        b = np.array([pvec[6], pvec[7], pvec[8]], dtype=float)  # linéaire
        c = float(pvec[9])                                      # constant
        return Ae, b, c

    def try_params(pvec):
        Ae, b, c = unpack(pvec)
        # Centre o = -0.5 * Ae^{-1} b
        try:
            o = -0.5 * np.linalg.solve(Ae, b)
        except np.linalg.LinAlgError:
            return None
        # Normalisation d'échelle: k tel que (x-o)^T (k Ae) (x-o) = 1
        denom = float(o.T @ Ae @ o - c)
        if abs(denom) < 1e-18:
            return None
        k = 1.0 / denom
        Ae_k = Ae * k
        # On exige Ae_k définie positive
        try:
            L = np.linalg.cholesky(Ae_k)   # L L^T = Ae_k
        except np.linalg.LinAlgError:
            return None
        A = L.T                            # A^T A = Ae_k
        return o, A, Ae_k

    out = try_params(p)
    if out is None:
        out = try_params(-p)
    if out is None:
        raise RuntimeError("Échec de l'ajustement de l'ellipsoïde.")

    o, A, Ae = out
    return o, A, Ae

def calibrate_xyz(raw_xyz: np.ndarray):
    o, A, Ae = fit_ellipsoid_ls(raw_xyz)
    # y = A @ (x - o) pour chaque échantillon
    Y = (raw_xyz - o) @ A.T
    return o, A, Y

def main(csv_path="datacalibration3.txt", out_path="mag_calibrated.csv"):
    csv_path = Path(csv_path)
    df = pd.read_csv(csv_path)

    # Utilise uniquement xmag, ymag, zmag
    M = df[["xmag", "ymag", "zmag"]].values.astype(float)

    # Calibrer
    o, A, Y = calibrate_xyz(M)

    # Résumés
    norms_before = np.linalg.norm(M, axis=1)
    norms_after = np.linalg.norm(Y, axis=1)

    print("Biais o (hard-iron):")
    print(o)  # shape (3,)
    print("\nMatrice de calibration A (y = A @ (x - o)) :")
    print(A)  # 3x3
    print("\nNorme moyenne avant:", norms_before.mean())
    print("Norme moyenne après :", norms_after.mean(),
          "(écart-type:", norms_after.std(), ")")

    # Sauvegarde des données calibrées
    out = df.copy()
    out["xmag_cal"] = Y[:, 0]
    out["ymag_cal"] = Y[:, 1]
    out["zmag_cal"] = Y[:, 2]
    out.to_csv(out_path, index=False)
    print(f"\nÉcrit: {out_path}")

if __name__ == "__main__":
    # Adapter les chemins si besoin
    main("datacalibration6.txt", "mag_calibrated.csv")
