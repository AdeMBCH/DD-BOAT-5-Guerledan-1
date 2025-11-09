# -*- coding: utf-8 -*-
# Script autonome de traçage proie vs bateau
# Python 3.5+

import math
import numpy as np
import matplotlib.pyplot as plt

# ---------- Paramètres ----------
Cx, Cy = -300.8, 59.9 #-177.5, 19.1                 # centre C [m]
R1 = 240                      # rayon proie [m]
w1 = -2.0*math.pi/(R1)   # vitesse angulaire proie [rad/s]
N = 12                          # nombre total de points
i = 5                             # index du point suivi -> phi0 = 2π*i/N
phi0 = 2.0*math.pi*(i/float(N))   # déphasage
T = 900                  # durée de simu [s]
dt = 0.1                          # pas [s]

# ---------- Fonctions ----------

def R2(t):
    # rayon d’orbite décroissant
    return 10.0*np.exp(-t/1000.0) + 5

def w2(t):
    return 2.0*math.pi/60.0

def prey_xy(t):
    # p(t) = C + R1 [cos(w1 t), sin(w1 t)]
    px = Cx + R1*np.cos(w1*t)
    py = Cy + R1*np.sin(w1*t)
    return px, py

def boat_xy(t):
    # a(t) = p(t) + R2(t) [cos(theta2), sin(theta2)]
    # theta2(t) = phi0 + w2(t)*t (modèle expérimental)
    px, py = prey_xy(t)
    r = R2(t)
    th = phi0 + w2(t)*t
    ax = px + r*np.cos(th)
    ay = py + r*np.sin(th)
    return ax, ay

# ---------- Simulation ----------
t = np.arange(0.0, T+dt, dt)
px, py = prey_xy(t)
ax, ay = boat_xy(t)

# ---------- Tracés ----------
plt.figure(figsize=(7,7))
plt.plot(px, py, label="Proie p(t)")
plt.plot(ax, ay, label="Bateau a(t)")
plt.scatter([px[0]],[py[0]], marker='o', zorder=3, label="Départ commun")
plt.axis('equal')
plt.xlabel("x Est [m]")
plt.ylabel("y Nord [m]")
plt.title("Trajectoires proie et bateau")
plt.grid(True)
plt.legend(loc="best")
plt.show()
