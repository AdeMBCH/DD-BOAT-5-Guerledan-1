"""Outil pedagogique pour déterminer et visualiser des trajectoires de proies et de prédateurs.
   Le code calcule les trajectoires dans un plan tangent local, génère des fichiers GPX pour les trajectoires
   et crée des visualisations statiques et animées des mouvements de la proie et des prédateurs.
   La partie affichage et géneration d'output .GPX sont réalisés avec l'aide de ChatGPT et peut demander des améliorations ou des adaptations.

   ##### Structure du code :#####
   Le code est structuré en sections ('régions') clairement délimitées pour faciliter la lecture.
   les sections sont :
    - Imports : importation des bibliothèques nécessaires
    - Fonctions de position : calcul des positions de la proie et des prédateurs
    - Fonctions et outils : fonctions auxiliaires pour la conversion GPX et l'animation
   Dans le main :
        - Paramètres : définition des paramètres de la simulation
        - Simulation principale : définition des paramètres, calcul des trajectoires et export GPX
        - Configuration des graphiques statiques 
        - Configuration des graphiques animés


##### Choix des paramètres :#####
   Les paramètres de la simulation sont définis au début du script. Ceux pertinants pour coller les contraintes physiques sont :
    ajustable selon le tracé voulu de la trajectoire de la proie (contrainre spatio-temporelle) :
        - durée : durée totale de la simulation en secondes
        - center_latlon : latitude et longitude du centre du plan tangent local en degrés
        - R1 : rayon de la trajectoire de la proie en mètres
            Sur le terrain il faut choisir R1 et center_latlon pour ne pas entrer en collision avec le rivage ou autre obstacle.
        
    ajustable en fonction de l'esthetique recherchée:
        - tau : constante de temps pour l'approche des prédateurs en secondes
            Pour que l'exponantielle converge, on voudra tau<durée/3. Ainsi on atteint 95% de la valeur finale R2i
        - R2i et R2e : Rayon interieur et exterieur désirée des prédateurs par rapport à la proie en mètres
            R2i<R2e; on prendra R2e en fonction de l'espace disponible physiquement et R2i pour permettre à tous les robots de ne pas se cogner.
    
    ajustable mais limitées en valeurs superieures par contraintes materielles :
        - speedmax_pred : vitesse maximale des prédateurs en m/s (limitée à 1.5 m/s)
            On tient toujours compte que la vitesse R2/omega2 est relative à la proie dont la vitesse est déjà speed_prey
        - speed_prey : vitesse de la proie en m/s 
            pour un bon résultat visuel, choisir moins de la moitié de speedmax_pred est un bon compromis.
        - N : nombre de prédateurs (selon les participants disponibles)
        - toff : décalage spatial exprimé comme un temps (en secondes) pour choisir le départ de la proie sur le cercle trajectoire
            ici utile pour éviter de commencer sous le pont

    déduits des paramètres précédents :
        - omega1 : vitesse angulaire de la proie en rad/s
        - omega2 : vitesse angulaire des prédateurs en rad/s autour de la proie

    A partir de ces paramètres, on peut visualiser clairement les trajectoires et justifier proprement des choix des paramètres omega. 
    Il est cependant necessaire de choisir 

    Le parametre t0 dans les fonctions de calcul des positions permet de tenir compte de la notion de temps de référence pour le départ.
    Inutile ici, il est necessaire pour une implémentation pour lequel on reglera l'heure de départ réelle.
"""

#region Imports
import math
import datetime
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import xml.etree.ElementTree as ET
#endregion


#region position functions

def prey(t, c, R1, omega1, t0, toff=0):
    """
    Position of the prey at a single time t. Returns a column vector (2, 1).
    """
    t = float(t)
    center = np.asarray(c, dtype=float).reshape(2, 1)
    angle = omega1 * (t - t0+toff)
    offset = np.array([[np.cos(angle)], [np.sin(angle)]])
    return center + R1 * offset

def target_pos(t, center_xy, R1, omega1, t0, tau, i, N):
    """
    Position of a predator at a single time t. Returns a column vector (2, 1).
    """
    base = prey(t, center_xy, R1, omega1, t0)               # computing from prey position
    R2 = (R2e-R2i) * np.exp(-(t - t0) / tau) + R2i
    omega2 = - (speedmax_pred-speed_prey)/R2
    phase_shift = 2 * np.pi * i / N                         # to get different positions for each predator
    angle = omega2 * (t - t0) + phase_shift
    offset = np.array([[np.cos(angle)], [np.sin(angle)]])
    return base + R2 * offset
#endregion

#region Functions and tools

#### Functions for GPX output ####
def xy_to_latlon(points_xy, center_latlon):
    """
    Convert local tangent plane (x east, y north in meters) to lat/lon.
    Used for GPX output.
    """
    lat0, lon0 = center_latlon
    latitudes = lat0 + points_xy[:, 1] / meters_per_deg_lat
    longitudes = lon0 + points_xy[:, 0] / meters_per_deg_lon
    return latitudes, longitudes

def write_prey_gpx(filename, latitudes, longitudes, times):
    """
    Write a GPX track for the prey trajectory.
    """
    gpx = ET.Element(
        "gpx",
        version="1.1",
        creator="tst.py",
        xmlns="http://www.topografix.com/GPX/1/1",
    )
    trk = ET.SubElement(gpx, "trk")
    ET.SubElement(trk, "name").text = "Prey Trajectory"
    trkseg = ET.SubElement(trk, "trkseg")
    for lat, lon, time_value in zip(latitudes, longitudes, times):
        trkpt = ET.SubElement(
            trkseg,
            "trkpt",
            lat=f"{lat:.8f}",
            lon=f"{lon:.8f}",
        )
        ET.SubElement(trkpt, "time").text = time_value.isoformat() + "Z"
    ET.ElementTree(gpx).write(filename, encoding="utf-8", xml_declaration=True)

#### Tools for animation ####
def _init():
    prey_line.set_data([], [])
    for line in target_lines:
        line.set_data([], [])
    return [prey_line, *target_lines]

def _update(frame_idx):
    idx = frame_idx
    prey_line.set_data(p[: idx + 1, 0], p[: idx + 1, 1])
    for line, M_traj in zip(target_lines, Ms):
        line.set_data(M_traj[: idx + 1, 0], M_traj[: idx + 1, 1])
    return [prey_line, *target_lines]

#endregion


if __name__ == "__main__":
    
    #region Parameters 
    durée = 750                                     # total duration in seconds
    center_latlon = (48.199706, -3.018784)          # latitude, longitude in degrees
    center_xy = np.array([[0.0], [0.0]])            # local tangent plane origin in meters
    speed_prey = 1/3                                # speed of prey in m/s
    speedmax_pred = 1/2+1/3                         # maximum speed of predators in m/s
    R1 = 240                                        # radius of prey trajectory 
    omega1 = -speed_prey/R1                          # angular speed of prey
    R2i = 5                                         # final distance of predators to prey 
    R2e = 15                                        # initial desired distance of predators to prey
    tau = 200                                       # time constant for predator approach
    omega2max = (speedmax_pred - speed_prey) / R2i  # maximum angular speed of predators in the prey reference frame [unusedin the code]
    N = 12                                          # number of predators   
    toff = 100                                      # time offset to avoid starting under the bridge
    # Conversion factors for lat/lon to local tangent plane for GPX output
    meters_per_deg_lat = 111320.0
    meters_per_deg_lon = 111320.0 * math.cos(math.radians(center_latlon[0]))
    #endregion
    
    ########

    #region Main simulation
    # showing needed parameters for ths simulation
    print (f'R1 = {R1}, R2 = {R2e}exp(-(t-t0)/{tau})+{R2i} omega1 = {speed_prey}/{R1}, omega1choisi=1/3R1 omega2max = {(speedmax_pred-speed_prey)}/R2 omega1choisi = 1/2R2')

    # compute trajectories and export GPX files if needed
    t_values = np.arange(0, durée, 0.1)
    fig_static, ax_static = plt.subplots()

    
    ###prey trajectory
    p = np.array([prey(ti, center_xy, R1, omega1, t0=0, toff=toff).flatten() for ti in t_values])
    
    ##GPX output for prey
    prey_latitudes, prey_longitudes = xy_to_latlon(p, center_latlon)
    start_time = datetime.datetime.utcnow().replace(microsecond=0)
    prey_times = [start_time + datetime.timedelta(seconds=float(seconds)) for seconds in t_values]
    write_prey_gpx("prey_trajectory.gpx", prey_latitudes, prey_longitudes, prey_times)

    ###predators trajectories
    Ms = []
    for idx in range(N):
        M = np.array([
                target_pos(ti, center_xy, R1, omega1, t0=0, tau=tau, i=idx,N=N,).flatten()
                for ti in t_values])
        Ms.append(M)
        ax_static.plot(M[:, 0], M[:, 1], label=f"M{idx} (target)")
        
        m_latitudes, m_longitudes = xy_to_latlon(M, center_latlon)
        start_time = datetime.datetime.utcnow().replace(microsecond=0)
        times = [
            start_time + datetime.timedelta(seconds=float(seconds))
            for seconds in t_values]
        write_prey_gpx(f"{idx+1}_trajectory.gpx", m_latitudes, m_longitudes, times)
    #endregion

    ############

    #region setting up the static plot
    ax_static.plot(p[:, 0], p[:, 1], color="black", lw=2, label="p (prey)")
    all_x = np.concatenate([p[:, 0]] + [traj[:, 0] for traj in Ms])
    all_y = np.concatenate([p[:, 1]] + [traj[:, 1] for traj in Ms])
    x_margin = 0.05 * (all_x.max() - all_x.min() or 1.0)
    y_margin = 0.05 * (all_y.max() - all_y.min() or 1.0)
    x_limits = (all_x.min() - x_margin, all_x.max() + x_margin)
    y_limits = (all_y.min() - y_margin, all_y.max() + y_margin)
    ax_static.set_xlabel("x position")
    ax_static.set_ylabel("y position")
    ax_static.set_title("Trajectories of p and M in the plane")
    ax_static.legend(ncol=3)
    ax_static.set_xlim(*x_limits)
    ax_static.set_ylim(*y_limits)
    ax_static.set_aspect("equal", adjustable="box")
    ax_static.grid(True)
    #endregion


    #region Setting up the animated plot
    fig_anim, ax_anim = plt.subplots()
    ax_anim.set_xlabel("x position")
    ax_anim.set_ylabel("y position")
    ax_anim.set_title("Animated Trajectories of p and M")
    ax_anim.grid(True)

    ax_anim.set_xlim(*x_limits)
    ax_anim.set_ylim(*y_limits)
    ax_anim.set_aspect("equal", adjustable="box")

    prey_line, = ax_anim.plot([], [], color="black", lw=2, label="p (prey)")
    target_lines = [ax_anim.plot([], [], label=f"M{idx} (target)")[0] for idx in range(N)]

    frame_indices = np.arange(len(t_values))
    frame_interval_ms = 100  # match the 0.1 s timestep

    anim = FuncAnimation(
        fig_anim,
        _update,
        frames=frame_indices,
        init_func=_init,
        blit=True,
        interval=frame_interval_ms,
        repeat=True,
    )
    ax_anim.legend(ncol=3)
    #endregion
    plt.show()
