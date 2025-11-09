from __future__ import print_function
import os
import sys
import time
import math
import argparse

import numpy as np

# Drivers
sys.path.append(os.path.join(os.path.dirname(__file__), 'drivers-ddboat-v2'))
import arduino_driver_v2 as arddrv
import imu9_driver_v2 as imudrv
import gps_driver_v2 as gpddrv

# Projections (compat Py3.5): Proj/transform
from pyproj import Proj
import gpxpy.gpx

# ---------- Paramètres repère / référence ----------

REF_LAT = 48.1991683
REF_LON = -3.01473

# REF_LAT = 48.19855
# REF_LON = -3.014635

# ---------- Outils ----------

def sur_360(a):
    return (a % 360.0 + 360.0) % 360.0

def ang_err(target, meas):
    # écart signé dans (-180, 180]
    return (target - meas + 180.0) % 360.0 - 180.0

# ---------- Conversion NMEA GLL DDMM.MMMM -> DD.DDDDD ----------

def gll_ddmm_to_dd(st):
    """st = [lat_ddmm.mmmm, N|S, lon_ddmm.mmmm, E|W, ...] -> (lat, lon) en degrés décimaux"""
    ilat = float(st[0])
    ns = st[1]
    ilon = float(st[2])
    ew = st[3]
    lat_deg = int(ilat / 100.0)
    lon_deg = int(ilon / 100.0)
    lat = lat_deg + (ilat - 100.0 * lat_deg) / 60.0
    lon = lon_deg + (ilon - 100.0 * lon_deg) / 60.0
    if ns == 'S':
        lat = -lat
    if ew == 'W':
        lon = -lon
    return lat, lon

# ---------- Calibration magnétomètre ----------
# Fichier: datacalibration.txt avec colonnes xmag, ymag, zmag

def _fit_ellipsoid_ls(M):
    x = M[:, 0]; y = M[:, 1]; z = M[:, 2]
    D = np.column_stack([x*x, y*y, z*z, x*y, x*z, y*z, x, y, z, np.ones_like(x)])
    U, S, Vt = np.linalg.svd(D, full_matrices=False)
    p = Vt[-1, :]

    def unpack(pv):
        Ae = np.array([
            [pv[0],     pv[3]/2.0, pv[4]/2.0],
            [pv[3]/2.0, pv[1],     pv[5]/2.0],
            [pv[4]/2.0, pv[5]/2.0, pv[2]]
        ], dtype=float)
        b = np.array([pv[6], pv[7], pv[8]], dtype=float)
        c = float(pv[9])
        return Ae, b, c

    def try_params(pv):
        Ae, b, c = unpack(pv)
        try:
            o = -0.5 * np.linalg.solve(Ae, b)
        except Exception:
            return None
        denom = float(np.dot(o, np.dot(Ae, o)) - c)
        if abs(denom) < 1e-18:
            return None
        k = 1.0 / denom
        Ae_k = Ae * k
        try:
            L = np.linalg.cholesky(Ae_k)
        except Exception:
            return None
        A = L.T
        return o, A

    out = try_params(p)
    if out is None:
        out = try_params(-p)
    if out is None:
        raise RuntimeError("Échec ajustement ellipsoïde")
    return out  # (o, A)


def load_mag_calibration(csv_path):
    if not os.path.isfile(csv_path):
        return None
    import csv
    with open(csv_path, newline='') as f:
        rdr = csv.DictReader(f)
        rows = list(rdr)
    if not rows:
        return None
    M = np.array([[float(r['xmag']), float(r['ymag']), float(r['zmag'])] for r in rows], dtype=float)
    o, A = _fit_ellipsoid_ls(M)
    return (A, o)

# ---------- Contrôleur de cap ----------

# class HeadingCtrl(object):
#     def __init__(self, imu, ard, A_mag=None, b_mag=None, kp=1.5, sat=230):
#         self.imu = imu
#         self.ard = ard
#         self.A_mag = A_mag
#         self.b_mag = b_mag
#         self.kp = kp
#         self.sat = sat

#     def read_heading_deg(self):
#         xraw, yraw, zraw = self.imu.read_mag_raw()
#         if (self.A_mag is not None) and (self.b_mag is not None):
#             v = np.array([xraw, yraw, zraw]) - self.b_mag
#             x, y, z = np.dot(self.A_mag, v)
#         else:
#             x, y = xraw, yraw
#         hdg = float(self.imu.heading_raw_deg(x, y))
#         return sur_360(hdg)

#     def step(self, hdg_cmd, vitesse):
#         hdg_meas = self.read_heading_deg()
#         print(hdg_meas)
#         err = ang_err(hdg_cmd, hdg_meas)
#         turn = self.kp * err
#         if turn > self.sat:
#             turn = self.sat
#         if turn < -self.sat:
#             turn = -self.sat
#         fwd = float(vitesse)
#         if fwd > self.sat:
#             fwd = self.sat
#         if fwd < -self.sat:
#             fwd = -self.sat
#         left  = int(max(-self.sat, min(self.sat, fwd + turn)))
#         right = int(max(-self.sat, min(self.sat, fwd - turn)))
#         self.ard.send_arduino_cmd_motor(left, right)
#         return err, hdg_meas, left, right
    
class HeadingCtrl(object):
    def __init__(self, imu, ard, A_mag=None, b_mag=None, kp=1.5, ki=0.5, kd=0, sat=230):
        self.imu, self.ard = imu, ard
        self.A_mag, self.b_mag = A_mag, b_mag
        self.kp, self.ki, self.kd, self.sat = kp, ki, kd, sat
        self.e_int = 0.0
        self.prev_err = 0.0
        self.prev_t = time.time()

        self.deadband_deg = 2.0         # bande morte cap
        self.v_cmd_prev = 0.0           # rampe vitesse
        self.prev_hdg = None            # pour D sur mesure
        self.turn_sign = +1             # +1 / -1 selon sens moteurs
        self.trim = 0.0                 # correction de biais propulseurs

    def read_heading_deg(self):
        xraw, yraw, zraw = self.imu.read_mag_raw()
        if (self.A_mag is not None) and (self.b_mag is not None):
            v = np.array([xraw, yraw, zraw]) - self.b_mag
            x, y, z = np.dot(self.A_mag, v)
        else:
            x, y = xraw, yraw
        hdg = float(self.imu.heading_raw_deg(x, y))
        return sur_360(hdg)

    # def step(self, hdg_cmd, vitesse):
    #     hdg_meas = self.read_heading_deg()
    #     print(hdg_meas)
    #     err = ang_err(hdg_cmd, hdg_meas) 
    #     t = time.time()
    #     dt = max(1e-2, t - self.prev_t)

    #     derr = (err - self.prev_err) / dt

    #     turn_unsat = self.kp*err + self.ki*self.e_int + self.kd*derr
    #     turn = max(-self.sat, min(self.sat, turn_unsat))
    #     if abs(turn) < self.sat*0.98:
    #         self.e_int += err * dt

    #     left  = int(max(-self.sat, min(self.sat, float(vitesse) + turn)))
    #     right = int(max(-self.sat, min(self.sat, float(vitesse) - turn)))
    #     self.ard.send_arduino_cmd_motor(left, right)

    #     self.prev_err, self.prev_t = err, t
    #     return err, hdg_meas, left, right

    def step(self, hdg_cmd, vitesse):
        hdg_meas = self.read_heading_deg()
        t = time.time()
        dt = max(1e-2, t - self.prev_t)

        # erreur cap avec bande morte
        err = ang_err(hdg_cmd, hdg_meas)
        if abs(err) < self.deadband_deg:
            err = 0.0

        # dérivée sur la mesure pour éviter le "kick"
        dmeas = 0.0 if self.prev_hdg is None else ang_err(hdg_meas, self.prev_hdg) / dt
        turn_unsat = self.turn_sign * (self.kp*err + self.ki*self.e_int - self.kd*dmeas)
        turn = max(-self.sat, min(self.sat, turn_unsat))

        # anti-windup simple
        if abs(turn) < self.sat*0.98:
            self.e_int += err * dt

        # gating vitesse selon l'erreur de cap
        E_TURN, E_SLOW = 30.0, 25.0
        gain = max(0.0, 1.0 - abs(err)/E_SLOW)
        v_eff = float(vitesse) * gain
        if abs(err) > E_TURN:
            v_eff = 0.0

        # rampe d'accélération
        slew = 50.0  # unités/s
        v_eff = max(min(v_eff, self.v_cmd_prev + slew*dt), self.v_cmd_prev - slew*dt)
        self.v_cmd_prev = v_eff

        # mixage priorisant le virage
        v_lim = max(0.0, self.sat - abs(turn))
        v_eff = max(min(v_eff, v_lim), -v_lim)

        left  = int(max(-self.sat, min(self.sat, v_eff + turn + self.trim)))
        right = int(max(-self.sat, min(self.sat, v_eff - turn - self.trim)))
        self.ard.send_arduino_cmd_motor(left, right)

        self.prev_t, self.prev_hdg = t, hdg_meas
        return err, hdg_meas, left, right

    
# ---------- Navigation point-à-point ----------

def bearing_geo_from_dxdy(dx, dy):
    # angle trigo [deg] depuis axe Est -> cap géographique 0=N, CW+
    ang_trigo_deg = math.degrees(math.atan2(dy, dx))
    hdg = 90.0 - ang_trigo_deg
    return sur_360(hdg)

class NavState(object):
    def __init__(self, x, y, lat, lon):
        self.x = x
        self.y = y
        self.lat = lat
        self.lon = lon

# Lecture GPS non bloquante avec NavState

def poll_nav(gps):
    ok, st = gps.read_gll_non_blocking()
    if not ok:
        return None
    lat, lon = gll_ddmm_to_dd(st)
    x, y = local_xy(lon,lat,REF_LON,REF_LAT)#wgs84_to_local(lon, lat)
    return NavState(x=x, y=y, lat=lat, lon=lon)

# def wait_first_fix(gps, timeout_s):
#     t0 = time.time()
#     last = None
#     while time.time() - t0 < timeout_s:
#         ns = poll_nav(gps)
#         if ns is not None:
#             return ns
#         time.sleep(0.2)
#     raise RuntimeError("GPS: pas de fix GLL dans le délai")

def wait_first_fix(gps, timeout_s):
    t0 = time.time()
    last = None
    while time.time() - t0 < timeout_s:
        ns = poll_nav(gps)
        if ns is not None:
            if last is None:
                last = ns
            else:
                dx = ns.x - last.x; dy = ns.y - last.y
                if abs(dx) + abs(dy) > 0.2:   # > 20 cm ⇒ fix frais
                    return ns
        time.sleep(0.2)
    raise RuntimeError("GPS: pas de fix GLL dans le délai")



# def navigate_to(ctrl, gps, goal_xy, vitesse_false, arrive_radius_m, timeout_s, log_seg):
#     x_goal, y_goal = goal_xy
#     t0 = time.time()
#     last_ns = wait_first_fix(gps, timeout_s=min(60.0, timeout_s))
#     while True:
#         ns = poll_nav(gps)
#         if ns is None:
#             ns = last_ns
#         else:
#             last_ns = ns
#             if log_seg is not None:
#                 log_seg.points.append(gpxpy.gpx.GPXTrackPoint(ns.lat, ns.lon))
#         dx = x_goal - ns.x
#         dy = y_goal - ns.y
#         dist = math.hypot(dx, dy)
#         if dist <= arrive_radius_m:
#             ctrl.ard.send_arduino_cmd_motor(0, 0)
#             return ns
#         hdg_cmd = bearing_geo_from_dxdy(dx, dy)
#         print(hdg_cmd)
#         ctrl.step(hdg_cmd, vitesse_false)
#         if time.time() - t0 > timeout_s:
#             ctrl.ard.send_arduino_cmd_motor(0, 0)
#             raise RuntimeError("Timeout navigation")
#         time.sleep(0.1)

def navigate_to(ctrl, gps, goal_xy, vitesse_false, arrive_radius_m, timeout_s, log_seg,lookahead_m=15.0):
    x_goal, y_goal = goal_xy
    t0 = time.time()
    # point A
    last_ns = wait_first_fix(gps, timeout_s=min(60.0, timeout_s))
    x_start, y_start = last_ns.x, last_ns.y

    # Géométrie du segment A->B
    dx_line = x_goal - x_start
    dy_line = y_goal - y_start
    line_len = math.hypot(dx_line, dy_line)
    if line_len < 1e-6:
        ctrl.ard.send_arduino_cmd_motor(0, 0)
        return last_ns
    vx = dx_line / line_len
    vy = dy_line / line_len

    while True:
        ns = poll_nav(gps)
        if ns is None:
            ns = last_ns
        else:
            last_ns = ns
            if log_seg is not None:
                log_seg.points.append(gpxpy.gpx.GPXTrackPoint(ns.lat, ns.lon))

        # Avancée le long de A->B et distance r à B
        dx = ns.x - x_start
        dy = ns.y - y_start
        s_along = dx * vx + dy * vy
        dist_goal = math.hypot(x_goal - ns.x, y_goal - ns.y)

        # Arrêt si atteint ou si on a dépassé la perpendiculaire de B
        if (dist_goal <= arrive_radius_m) or (s_along >= line_len):
            ctrl.ard.send_arduino_cmd_motor(0, 0)
            return ns

        # Calcul cap
        # Point de visée L mètres devant la projection
        s_ref = s_along + float(lookahead_m)
        if s_ref < 0.0:
            s_ref = 0.0
        if s_ref > line_len:
            s_ref = line_len
        x_ref = x_start + vx * s_ref
        y_ref = y_start + vy * s_ref
        hdg_cmd = bearing_geo_from_dxdy(x_ref - ns.x, y_ref - ns.y)
        print(hdg_cmd)
        ctrl.step(hdg_cmd, vitesse_false)

        if time.time() - t0 > timeout_s:
                ctrl.ard.send_arduino_cmd_motor(0, 0)
                raise RuntimeError("Timeout navigation")
        time.sleep(0.1)

def local_xy(long_deg, lat_deg, long0_deg, lat0_deg, R=6378137.0):
    long = math.radians(long_deg)
    lat = math.radians(lat_deg)
    long0 = math.radians(long0_deg)
    lat0 = math.radians(lat0_deg)

    dlong = (long - long0)
    dlat = lat - lat0
    x = R * math.cos(lat0) * dlong
    y = R * dlat
    return x, y


# ---------- Main ----------

def main():
    p = argparse.ArgumentParser(description="Aller-retour vers un waypoint depuis Guérlédan")
    grp = p.add_mutually_exclusive_group(required=True)
    grp.add_argument('--lat', type=float, help='Latitude cible en degrés décimaux')
    p.add_argument('--lon', type=float, help='Longitude cible en degrés décimaux')
    grp.add_argument('--xy', nargs=2, type=float, metavar=('X', 'Y'), help='Cible locale en mètres [Est, Nord]')
    p.add_argument('--vitesse_false', type=float, default=60.0, help='Consigne vitesse [-100,100]')
    p.add_argument('--arrive_radius', type=float, default=2.0, help="Rayon d'arrivée en mètres")
    p.add_argument('--timeout', type=float, default=900.0, help='Timeout total par jambe [s]')
    p.add_argument('--gpx', type=str, default='trace_nav.gpx', help='Fichier GPX de trace')
    args = p.parse_args()

    if (args.lat is not None) and (args.lon is None):
        p.error('--lat nécessite --lon')

    # Cible locale
    if args.xy is not None:
        goal_xy = (float(args.xy[0]), float(args.xy[1]))
    else:
        goal_xy = local_xy(args.lon, args.lat,REF_LON,REF_LAT) #wgs84_to_local(args.lon, args.lat)

    # Périphériques
    ard = arddrv.ArduinoIO()
    imu = imudrv.Imu9IO()
    gps = gpddrv.GpsIO(tty_dev=0)
    try:
        gps.set_filter_speed("0")
    except Exception:
        pass

    # Calibration magnétomètre
    calib = load_mag_calibration("datacalibration.txt")
    A_mag, b_mag = calib

    ctrl = HeadingCtrl(imu=imu, ard=ard, A_mag=A_mag, b_mag=b_mag)

    # GPX log
    gpx = gpxpy.gpx.GPX()
    trk = gpxpy.gpx.GPXTrack()
    gpx.tracks.append(trk)
    seg = gpxpy.gpx.GPXTrackSegment()
    trk.segments.append(seg)

    # HOME
    # home = wait_first_fix(gps, timeout_s=60.0)
    # print("HOME lat=%.6f lon=%.6f x=%.1f y=%.1f" % (home.lat, home.lon, home.x, home.y))

    # HOME = origine Guérlédan (repère fixe)
    home = NavState(x=0.0, y=0.0, lat=REF_LAT, lon=REF_LON)
    print("HOME lat=%.6f lon=%.6f x=0.0 y=0.0" % (home.lat, home.lon))

    # ALLER
    print("Aller vers x=%.1f y=%.1f vitesse=%.1f" % (goal_xy[0], goal_xy[1], args.vitesse_false))
    last_at_goal = navigate_to(ctrl, gps, goal_xy, args.vitesse_false, args.arrive_radius, args.timeout, seg)
    print("Arrivé proche cible: lat=%.6f lon=%.6f" % (last_at_goal.lat, last_at_goal.lon))

    # Pause
    ctrl.ard.send_arduino_cmd_motor(0, 0)
    time.sleep(2.0)

    # RETOUR HOME
    print("Retour HOME")
    _ = navigate_to(ctrl, gps, (home.x, home.y), args.vitesse_false, args.arrive_radius, args.timeout, seg)

    ctrl.ard.send_arduino_cmd_motor(0, 0)

    # Sauvegarde GPX
    with open(args.gpx, 'w') as f:
        f.write(gpx.to_xml())
    print("Trace GPX -> %s" % args.gpx)

if __name__ == '__main__':
    try:
        main()
    finally:
        # fail-safe stop en toute circonstance
        try:
            arddrv.ArduinoIO().send_arduino_cmd_motor(0, 0)
        except Exception:
            pass
