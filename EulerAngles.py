import csv
import os
import sys
import time
from pathlib import Path
from math import pi
import numpy as np


DATA_FILE = "datacalibration.txt"
G0 = 9.81  # for normalizing accelerometer data

#region Ellipsoid fitting and calibration

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
#endregion

#region Main calibration function 
def Calibration(csv_path="datacalibration.txt", out_path="mag_calibrated.csv"):
   csv_path = Path(csv_path)
   out_path = Path(out_path)

   with csv_path.open(newline="") as f:
       reader = csv.DictReader(f)
       fieldnames = reader.fieldnames or []
       rows = list(reader)

   # Utilise uniquement xmag, ymag, zmag
   M = np.array(
       [[float(row["xmag"]), float(row["ymag"]), float(row["zmag"])] for row in rows],
       dtype=float,
   )

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
   extra_fields = ["xmag_cal", "ymag_cal", "zmag_cal"]
   out_fieldnames = fieldnames + [name for name in extra_fields if name not in fieldnames]

   with out_path.open("w", newline="") as f:
       writer = csv.DictWriter(f, fieldnames=out_fieldnames)
       writer.writeheader()
       writer.writerows(rows)

   return A, o
#endregion

#region tools

def Cal(mag, A_mag, b_mag):
    x = (np.asarray(mag).reshape(3,)- b_mag)
    return A_mag @ x

def ang_diff(a, b):
    # écart signé dans (-180, 180]
    return (a - b + 180.0) % 360.0 - 180.0
#endregion

#region euler angles from accelerometer and magnetometer

def R(a, y):
    a = np.array(a)/np.linalg.norm(a)
    y = np.array(y)
    n=(y - (y.T@a)*a)
    d = np.linalg.norm(n)
    return np.array([n/d, np.cross(a,n/d), a])


def eulermat2angles(R):
    φ=np.arctan2(R[2,1],R[2,2])
    θ=-np.arcsin(R[2,0])
    ψ=np.arctan2(R[1,0],R[0,0])
    return φ,θ,ψ

#test with a boat flat and heading north

a0=np.array([[0],[0],[1]])
I=1.15

y0=np.array([[np.cos(I)],[0],[-np.sin(I)]])

print(R(a0,y0))

#endregion


#region extract data and prepare mission

A_mag, b_mag = Calibration("datacalibration.txt", "mag_calibrated.csv")
sys.path.append(os.path.join(os.path.dirname(__file__), 'drivers-ddboat-v2'))
import arduino_driver_v2 as arddrv
import imu9_driver_v2 as imudrv
ard = arddrv.ArduinoIO()
imu = imudrv.Imu9IO()

xmag,ymag,zmag = imu.read_mag_raw()
mag = [[xmag,ymag,zmag]]
A_mag, b_mag = main()
#np.set_printoptions(precision=6, suppress=True)
print("A_mag:\n", A_mag)
print("bneg_mag:", b_mag)
#endregion

#region mission
### Route to 270° ### this mission is an old version for heading control and is now replaced. It however shows a functional example 
                    # of using the calibration and the Euler angles adjusted with accelerometers.
                    # While being coherent when motors are off, is shows a lot of noise when motors are running due to vibrations affecting the accelerometers.

heading = 270.0

t0 = time.time()


while time.time() - t0 < 180.0:
    
    
    xmagr, ymagr, zmagr = imu.read_mag_raw()            # raw magnetometer readings
    X = np.array([[xmagr],[ymagr],[zmagr]])             # arranged in a column vector
    y1 = (Cal(X, A_mag, b_mag))                         # calibrated magnetometer readings    
    a1 = imu.read_accel_raw()                           # raw accelerometer readings    
    
    phi, theta, psi = eulermat2angles(R(a1,y1))         # Euler angles from accelerometer and magnetometer
    heading_meas = imu.heading_raw_deg(y1[0], y1[1])    # heading from magnetometer only
    
    phi, theta, psi = phi*180.0/pi, theta*180.0/pi, (psi*180.0/pi) + 180
                                                        # adjust psi to [0,360]
    print (phi, theta, psi, heading_meas)               # display Euler angles and heading for comparison
    

    # proportional controller for heading
    # take into account the dead zone
    # still show some defects on the chosen values >>> we need to get a feed forward control to handle the imprecision
    # a PID will help with the correct input needed for FFC
    

    delta = ang_diff(heading, psi)
    absdelta = abs(delta)

    if absdelta > 135.0 :
        acc = 100 + (absdelta-135)*100/45 
        dec = -acc
    elif absdelta > 90.0 :
        acc = 200 - (absdelta)*100/135
        dec = -100
    elif absdelta > 0.0:
        acc = 200 - (absdelta)*100/135
        dec = 200 - (absdelta)*150/90
    else:
        acc = +100
        dec = -100


    if delta < 0.0 :          #aller à droite
        left_speed = acc
        right_speed = dec
    else:                    #aller à gauche
        left_speed = dec
        right_speed = acc


    #ard.send_arduino_cmd_motor(left_speed,right_speed) ### uncomment to enable motor commands
    time.sleep(0.1)
ard.send_arduino_cmd_motor(0,0) # stop
#endregion
