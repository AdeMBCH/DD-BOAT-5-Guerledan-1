import os
import sys
import time
import numpy as np

DATA_FILE = "datacalibration.txt"
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

def RMS(D):
    D = D.astype(float, copy=False)
    b = np.ones(D.shape[0], dtype=float)
    return np.linalg.lstsq(D, b, rcond=-1.0)[0]

def inv_(Q):
  w, V = np.linalg.eigh(0.5*(Q+Q.T))
  tol = max(1e-12, 1e-6*np.max(np.abs(w)))
  w_clip = np.where(w > tol, w, tol)
  Qinv = V @ np.diag(1.0 / w_clip) @ V.T
  return Qinv

def A_et_b(X, scale=1.0):
  p = RMS(M(X))
  p1,p2,p3,p4,p5,p6,p7,p8,p9,p10 = p
  Q = np.array([[p1, 0.5*p4, 0.5*p5],
                [0.5*p4, p2,  0.5*p6],
                [0.5*p5, 0.5*p6, p3]], float)
  ell = np.array([p7, p8, p9], float)
  Q = 0.5*(Q+Q.T)
  
  Qinv = inv_(Q)
  b = -0.5 * (Qinv @ ell)

 #normalisation
  k = 1.0 - p10 + 0.25 * (ell @ (Qinv @ ell))
  if k <= 0.0:
      # corrige le signe global si nécessaire
      Q  = -Q
      ell = -ell
      p10 = -p10
      Qinv = inv_(Q)
      b = -0.5 * (Qinv @ ell)
      k = 1.0 - p10 + 0.25 * (ell @ (Qinv @ ell))

  Q_unit = Q / k

  # A = scale * sqrt(Q_unit)
  w, V = np.linalg.eigh(Q_unit)
  w = np.maximum(w, 1e-12)
  A = scale * (V @ np.diag(np.sqrt(w)) @ V.T)

  return A, -b

mag = load_mag()
A_mag, bneg_mag = A_et_b(mag, scale=1.0)
print("A_mag:\n", A_mag)
print("bneg_mag:", bneg_mag)

####

sys.path.append(os.path.join(os.path.dirname(__file__), 'drivers-ddboat-v2'))
import arduino_driver_v2 as arddrv
import imu9_driver_v2 as imudrv
ard = arddrv.ArduinoIO()
imu = imudrv.Imu9IO()


def Cal(mag, A_mag, bneg_mag):
    x = (np.asarray(mag).reshape(3,) + bneg_mag)
    return A_mag @ x  


### Route to 270°


Kpa = 10
Kpd = 100
heading = 270.0


t0 = time.time()
def ang_diff(a, b):
    # écart signé dans (-180, 180]
    return (a - b + 180.0) % 360.0 - 180.0

heading = 270.0
while time.time() - t0 < 50.0:
    xm, ym, zm = imu.read_mag_raw()
    xmag, ymag, zmag = Cal([xm, ym, zm], A_mag, bneg_mag)
    heading_meas = imu.heading_raw_deg(xmag, ymag)

    delta = ang_diff(heading, heading_meas)

    acc, dec = 100, 0
    if delta > 0:          # besoin de tourner à GAUCHE (CCW)
        left_speed  = dec
        right_speed = acc
    else:                   # besoin de tourner à DROITE (CW)
        left_speed  = acc
        right_speed = dec

    ard.send_arduino_cmd_motor(left_speed, right_speed)
    time.sleep(0.01)

ard.send_arduino_cmd_motor(0, 0)

