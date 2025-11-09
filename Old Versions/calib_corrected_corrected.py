# calibrate_magnetometer_nopandas.py
# Ajustement ellipsoïde LS -> biais o et matrice de gain A (3D), puis y = A @ (x - o)
# Entrée: CSV avec en-tête contenant xmag, ymag, zmag
# Sortie: CSV identique + colonnes xmag_cal, ymag_cal, zmag_cal

import argparse
import csv
from pathlib import Path
import numpy as np

def fit_ellipsoid_ls(X: np.ndarray):
    """
    X: Nx3 array des mesures brutes [x, y, z]
    Retourne:
      o  : vecteur 3x1 biais (hard-iron)
      A  : matrice 3x3 telle que y = A @ (x - o) et ||y|| ≈ 1
      Ae : 3x3 symétrique avec A^T A = Ae
    """
    x, y, z = X[:, 0], X[:, 1], X[:, 2]

    # Matrice de design D avec d(x) = [x^2, y^2, z^2, x*y, x*z, y*z, x, y, z, 1]
    D = np.column_stack([x*x, y*y, z*z, x*y, x*z, y*z, x, y, z, np.ones_like(x)])

    # Résoudre D p = 0 par SVD (p = vecteur singulier droit associé à la plus petite SV)
    _, _, Vt = np.linalg.svd(D, full_matrices=False)
    p = Vt[-1, :]  # 10 paramètres

    def unpack(pvec):
        Ae = np.array([
            [pvec[0],       pvec[3]/2.0, pvec[4]/2.0],
            [pvec[3]/2.0,   pvec[1],     pvec[5]/2.0],
            [pvec[4]/2.0,   pvec[5]/2.0, pvec[2]]
        ], dtype=float)
        b = np.array([pvec[6], pvec[7], pvec[8]], dtype=float)
        c = float(pvec[9])
        return Ae, b, c

    def try_params(pvec):
        Ae, b, c = unpack(pvec)
        # Centre o = -0.5 * Ae^{-1} b
        try:
            o = -0.5 * np.linalg.solve(Ae, b)
        except np.linalg.LinAlgError:
            return None
        # Échelle k telle que (x-o)^T (k Ae) (x-o) = 1
        denom = float(o.T @ Ae @ o - c)
        if abs(denom) < 1e-18:
            return None
        k = 1.0 / denom
        Ae_k = Ae * k
        # Exiger Ae_k définie positive
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
    Y = (raw_xyz - o) @ A.T
    return o, A, Y

def read_csv_no_pandas(path: Path):
    """
    Lit un CSV avec en-tête. Retourne:
      header: liste des noms de colonnes
      rows  : liste de lignes (chaînes)
    """
    with path.open("r", newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        rows = [row for row in reader if row]  # ignorer lignes vides
    return header, rows

def write_csv_no_pandas(path: Path, header, rows):
    with path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)

def extract_mag_matrix(header, rows):
    try:
        ix = header.index("xmag")
        iy = header.index("ymag")
        iz = header.index("zmag")
    except ValueError:
        raise RuntimeError("Colonnes xmag, ymag, zmag introuvables dans l'en-tête.")
    M = np.empty((len(rows), 3), dtype=float)
    for i, r in enumerate(rows):
        M[i, 0] = float(r[ix])
        M[i, 1] = float(r[iy])
        M[i, 2] = float(r[iz])
    return M

def main(in_path="datacalibration.txt", out_path="mag_calibrated.csv"):
    in_path = Path(in_path)
    out_path = Path(out_path)

    # Lecture brute
    header, rows = read_csv_no_pandas(in_path)

    # Extraction magnétomètre
    M = extract_mag_matrix(header, rows)

    # Calibrage
    o, A, Y = calibrate_xyz(M)

    # Résumés
    norms_before = np.linalg.norm(M, axis=1)
    norms_after  = np.linalg.norm(Y, axis=1)

    print("Biais o (hard-iron):")
    print(o)
    print("\nMatrice de calibration A (y = A @ (x - o)) :")
    print(A)
    print("\nNorme moyenne avant:", float(norms_before.mean()))
    print("Norme moyenne après :", float(norms_after.mean()),
          "(écart-type:", float(norms_after.std()), ")")

    # Écriture CSV de sortie: colonnes d'origine + calibrées
    new_header = list(header) + ["xmag_cal", "ymag_cal", "zmag_cal"]
    out_rows = []
    for i, r in enumerate(rows):
        out_rows.append(list(r) + [f"{Y[i,0]:.9f}", f"{Y[i,1]:.9f}", f"{Y[i,2]:.9f}"])
    write_csv_no_pandas(out_path, new_header, out_rows)
    print(f"\nÉcrit: {out_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calibration magnétomètre sans pandas.")
    parser.add_argument("-i", "--input", default="datacalibration.txt",
                        help="Chemin du fichier CSV d'entrée")
    parser.add_argument("-o", "--output", default="mag_calibrated.csv",
                        help="Chemin du CSV de sortie")
    args = parser.parse_args()
    main(args.input, args.output)
