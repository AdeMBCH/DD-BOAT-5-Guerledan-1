# Autonomous_proie.py
from __future__ import print_function
import os
import sys
import time
import math
import argparse
import numpy as np
import signal

# Drivers
sys.path.append(os.path.join(os.path.dirname(__file__), 'drivers-ddboat-v2'))
import arduino_driver_v2 as arddrv
import imu9_driver_v2 as imudrv
import gps_driver_v2 as gpddrv

# Projections / GPX
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
    return (target - meas + 180.0) % 360.0 - 180.0  # [-180,180]

PI = math.pi
TWOPI = 2.0 * math.pi

def wrap_to_pi(a): 
    return (a + PI) % TWOPI - PI  # (-pi, pi]

def deg_to_rad(d): 
    return d * PI / 180.0

def rad_to_deg(r): 
    return r * 180.0 / PI

def bearing_geo_from_vec(vx, vy):
    # 0° = Nord, CW+. vx=Est, vy=Nord
    return (90.0 - math.degrees(math.atan2(vy, vx))) % 360.0

def wait_for(t0):
   """
   return the time.time() value when t0 is attained


   t0 is a str : "hh:mm:ss"
   """
   if not(type(t0)==str and len(t0)==8):
       print ('wrong type')
       return False
   if (not t0[0:2].isdigit()) or (not t0[3:5].isdigit()) or (not t0[6:8].isdigit()):
       print('wrong format')
       return False
   h0=int(t0[0:2])
   m0=int(t0[3:5])
   s0=int(t0[6:8])
   if not(0<=h0<24 and 0<=m0<60 and 0<=s0<60):
       print('not rights time numbers')
       return False
   t=time.strftime('%X')
   h=int(t[0:2])
   m=int(t[3:5])
   s=int(t[6:8])
   print("start in " + str(h - h0) + "h" + str(m - m0) + "m" + str(s - s0) + "s")
   while t != t0 and h*3600+m*60+s<=h0*3600+m0*60+s0 :
       t=time.strftime('%X')
       h=int(t[0:2])
       m=int(t[3:5])
       s=int(t[6:8])
       dt=-(h*3600+m*60+s)+(h0*3600+m0*60+s0)
       dh = dt//3600
       dm = dt%3600//60
       ds = dt%3600%60
       sys.stdout.write("\r" + " " * 20 + "\r")
       sys.stdout.write("start in " + str(dh) + "h" + str(dm) + "m" + str(ds) + "s")
       sys.stdout.flush()
       time.sleep(0.1)


   return time.time()


# ---------- Conversion NMEA ----------
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

# ---------- Calibration magnétomètre (CSV: xmag, ymag, zmag) ----------
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

    out = try_params(p) or try_params(-p)
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

# ---------- FFC ----------
class FFCtrl(object):
    def __init__(self, imu, ard, A_mag=None, b_mag=None, alpha_v=1.0, alpha_w=0.7, kp_psi=1.2, mmax=230, lam=3.0):
        self.imu = imu
        self.ard = ard

        self.A_mag = A_mag
        self.b_mag = b_mag

        self.alpha_v = float(alpha_v)
        self.alpha_w = float(alpha_w)

        self.kp_psi = float(kp_psi)
        self.mmax = int(mmax)
        self.lam = float(lam)

    def tanh_sat(self, x):# Sortie entre 200 et -200 et pas entre 1 et -1 pour matcher avec les moteurs
        return self.mmax * math.tanh(x / self.lam)

    def read_heading_rad(self):
        # Lecture brute + calibration ellipsoïde + cap géographique
        xraw, yraw, zraw = self.imu.read_mag_raw()
        if (self.A_mag is not None) and (self.b_mag is not None):
            v = np.array([xraw, yraw, zraw]) - self.b_mag
            x, y, _ = np.dot(self.A_mag, v)
        else:
            x, y = xraw, yraw
        hdg_deg = float(self.imu.heading_raw_deg(x, y))
        print(hdg_deg)
        return deg_to_rad(sur_360(hdg_deg))

    def step(self, u_d, omega_d, psi_d):
        # FFC pur avec correction P en cap, saturation tanh
        psi = self.read_heading_rad()
        e_psi = wrap_to_pi(psi - psi_d) # sawtooth
        omega_cmd = omega_d - self.kp_psi * e_psi

        left  = int(round(self.tanh_sat(self.alpha_v*u_d + self.alpha_w*omega_cmd)))
        right = int(round(self.tanh_sat(self.alpha_v*u_d - self.alpha_w*omega_cmd)))

        left = max(-self.mmax, min(self.mmax, left))
        right = max(-self.mmax, min(self.mmax, right))

        self.ard.send_arduino_cmd_motor(left, right)

        return left, right, e_psi, omega_cmd

# ---------- Navigation ----------
class NavState(object):
    def __init__(self, x, y, lat, lon):
        self.x = x
        self.y = y
        self.lat = lat
        self.lon = lon

# Lecture GPS non bloquante
def bearing_geo_from_dxdy(dx, dy):
    # angle trigo [deg] depuis axe Est -> cap géographique 0=N, CW+
    ang_trigo_deg = math.degrees(math.atan2(dy, dx))
    hdg = 90.0 - ang_trigo_deg
    return sur_360(hdg)

def poll_nav(gps):
    ok, st = gps.read_gll_non_blocking()
    if not ok:
        return None
    lat, lon = gll_ddmm_to_dd(st)
    x, y = local_xy(lon, lat, REF_LON, REF_LAT)
    return NavState(x=x, y=y, lat=lat, lon=lon)

def wait_first_fix(gps, timeout_s):
    t0 = time.time(); last = None
    while time.time() - t0 < timeout_s:
        ns = poll_nav(gps)
        if ns is not None:
            if last is None:
                last = ns
            else:
                dx = ns.x - last.x; 
                dy = ns.y - last.y
                if abs(dx) + abs(dy) > 0.2: # > 20 cm
                    last =ns  
                else:
                    return ns
        time.sleep(0.2)
    raise RuntimeError("GPS: pas de fix GLL dans le délai")

def goto_point(ctrl, gps, goal_xy, u_d, radius_m, timeout_s, seg=None):
    xg, yg = goal_xy
    t0 = time.time()
    last = wait_first_fix(gps, timeout_s=min(30.0, timeout_s))
    while True:
        ns = poll_nav(gps)
        if ns is None:
            ns = last
        else:
            last = ns
            if seg is not None:
                seg.points.append(gpxpy.gpx.GPXTrackPoint(ns.lat, ns.lon))

        dx, dy = xg - ns.x, yg - ns.y
        dist = math.hypot(dx, dy)
        if dist <= radius_m:
            ctrl.ard.send_arduino_cmd_motor(0, 0)
            return ns
        print(bearing_geo_from_vec(dx, dy))
        psi_d = deg_to_rad(bearing_geo_from_vec(dx, dy))
        ctrl.step(u_d=u_d, omega_d=0.0, psi_d=psi_d)

        if time.time() - t0 > timeout_s:
            ctrl.ard.send_arduino_cmd_motor(0, 0)
            raise RuntimeError("Timeout goto_point")
        time.sleep(0.1)

# pour rétrocompatibilité
def go_to_point(ctrl, gps, goal_xy, u_d, radius_m, timeout_s, seg=None):
    return goto_point(ctrl, gps, goal_xy, u_d, radius_m, timeout_s, seg)

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

# ---------- Trajectoire proie + orbite ----------
class PreyOrbit(object):
    def __init__(self, Cx, Cy, R1, w1, N, kidx=5):
        self.Cx = float (Cx)
        self.Cy = float (Cy)
        self.R1 = float(R1)
        self.w1 = float(w1)
        self.N = int(N)
        self.kidx = int(kidx)
        self.phi0 = 2.0*math.pi * (self.kidx / float(self.N))

    @staticmethod
    def R2(t):
        return 10.0*np.exp(-t/200.0) + 5

    def w2(self, t):
        return 1.0/(2*self.R2(t))
    
    def prey_pos(self, t):
        return (self.Cx + self.R1 * math.cos(self.w1 * (t+150)), self.Cy + self.R1 * math.sin(self.w1 * (t+150)))

    def desired_pos(self, t):
        px, py = self.prey_pos(t)
        R2t = self.R2(t)
        w2 = self.w2(t)
        return (px,py)
        # return (px + R2t * math.cos(w2*t+ 2.0 *math.pi*5/self.N),
        #         py + R2t * math.sin(w2*t+ 2.0 *math.pi*5/self.N))

    def desired_vel(self, t, dt=0.1):
        x1, y1 = self.desired_pos(t - dt)
        x2, y2 = self.desired_pos(t + dt)
        vx = (x2 - x1) / (2.0 * dt)
        vy = (y2 - y1) / (2.0 * dt)
        return vx, vy

# ---------- Suivi de a(t) ----------
def track_at(ctrl, gps, traj, t_start, u_d_nom, k_pos=0.8, dt_cmd=0.1, duration_s=750.0, seg_robot=None, seg_ref=None):
    
    # v_ref = dot(a,t) + k_pos*(a - x)

    last = wait_first_fix(gps, timeout_s=60.0)
    if seg_robot is not None:
        seg_robot.points.append(gpxpy.gpx.GPXTrackPoint(last.lat, last.lon))

    t0 = t_start
    psi_d_prev = None
    t_prev = None
    t_begin = time.time()

    while True:
        ns = poll_nav(gps)
        if ns is None:
            ns = last
        else:
            last = ns
            if seg_robot is not None:
                seg_robot.points.append(gpxpy.gpx.GPXTrackPoint(ns.lat, ns.lon))

        t = time.time() - t0
        xd, yd = traj.desired_pos(t)
        vxd, vyd = traj.desired_vel(t)

        ex, ey = xd - ns.x, yd - ns.y
        vrefx = vxd + k_pos * ex
        vrefy = vyd + k_pos * ey

        psi_d = deg_to_rad(bearing_geo_from_vec(vrefx, vrefy))
        if psi_d_prev is None:
            omega_d = 0.0
        else:
            dpsi = wrap_to_pi(psi_d - psi_d_prev)
            dt = max(1e-3, (time.time() - t_prev))
            omega_d = dpsi / dt

        ctrl.step(u_d=u_d_nom, omega_d=omega_d, psi_d=psi_d)

        psi_d_prev = psi_d
        t_prev = time.time()

        if time.time() - t_begin > duration_s:
            ctrl.ard.send_arduino_cmd_motor(0, 0)
            return
        time.sleep(max(0.02, dt_cmd))

# ---------- Main ----------
def main():
    duree=750.0 #trois endroits à changer

    p = argparse.ArgumentParser(description="Suivi proie connue avec FFC (tanh + sawtooth)")

    # Proie + orbite
    p.add_argument('--C_latlon', nargs=2, type=float, default=[48.199706,-3.018784],help='Centre C de la proie (deg)')
    p.add_argument('--R1', type=float, default=240, help='Rayon cercle de la proie [m]')
    p.add_argument('--w1', type=float, default=-1.0/(3*240), help='Vitesse angulaire proie [rad/s]')
    p.add_argument('--N', type=int, default= 9, help='N pour le déphasage 2π*5/N')
    p.add_argument('--kidx', type=int, default=5, help='Index du bateau')

    # FFC / guidance
    p.add_argument('--u_d', type=float, default=2, help='Avance normalisée [-1,1]')
    p.add_argument('--alpha_v', type=float, default=1.0)
    p.add_argument('--alpha_w', type=float, default=0.7)
    p.add_argument('--kp_psi', type=float, default=1.2)
    p.add_argument('--lam', type=float, default=3.0)
    p.add_argument('--mmax', type=int, default=200)
    p.add_argument('--k_pos', type=float, default=0.8, help='Gain position du champ de guidage')
    p.add_argument('--dt_cmd', type=float, default=0.1, help='Période de commande [s]')
    p.add_argument('--duration', type=float, default=duree, help='Durée de suivi [s]')

    p.add_argument('--arrive_radius', type=float, default=3.0)
    p.add_argument('--timeout', type=float, default=duree)
    p.add_argument('--gpx', type=str, default='trace_nav.gpx')

    # Attente de départ
    p.add_argument('--start_time', type=str, default="00:00:00", help="attente de  l'horaire 'hh:mm:ss")
    p.add_argument('--wait_key', action='store_true', help="Attendre Enter avant le suivi")
    p.add_argument('--start_after', type=float, default=0.0, help="Attente (s) avant suivi")

    args = p.parse_args()

    # IO
    ard = arddrv.ArduinoIO()
    imu = imudrv.Imu9IO()
    gps = gpddrv.GpsIO(tty_dev=0)
    try:
        gps.set_filter_speed("0")
    except Exception:
        pass

    # Calibration magnétomètre
    calib = load_mag_calibration("datacalibration.txt")
    if calib is None:
        A_mag = None; b_mag = None
    else:
        A_mag, b_mag = calib

    # Init FFC
    ctrl = FFCtrl(imu=imu, ard=ard, A_mag=A_mag, b_mag=b_mag,
                  alpha_v=args.alpha_v, alpha_w=args.alpha_w,
                  kp_psi=args.kp_psi, mmax=args.mmax, lam=args.lam)

    # GPX
    gpx = gpxpy.gpx.GPX()
    trk_robot = gpxpy.gpx.GPXTrack(); gpx.tracks.append(trk_robot)
    seg_robot = gpxpy.gpx.GPXTrackSegment(); trk_robot.segments.append(seg_robot)

    # Paramètres proie
    Cx, Cy = local_xy(args.C_latlon[1], args.C_latlon[0], REF_LON, REF_LAT)
    traj = PreyOrbit(Cx=Cx, Cy=Cy, R1=args.R1, w1=args.w1, N=args.N, kidx=args.kidx)

    # 1) HOME (fix GPS)
    # home = wait_first_fix(gps, timeout_s=60.0)
    # print("HOME lat=%.6f lon=%.6f" % (home.lat, home.lon))
    home = NavState(x=0.0, y=0.0, lat=REF_LAT, lon=REF_LON)
    print("HOME lat=%.6f lon=%.6f x=0.0 y=0.0" % (home.lat, home.lon))

    # # 2) Aller au point de départ de la proie p(0) = C + [R1, 0]
    # p0x, p0y = traj.prey_pos(0.0)
    # print("Goto p(0) -> x=%.1f y=%.1f" % (p0x, p0y))
    # _ = goto_point(ctrl, gps, (p0x, p0y), u_d=args.u_d,
    #                radius_m=args.arrive_radius, timeout_s=args.timeout, seg=seg_robot)

    # # 3) Attente de départ puis suivi
    # ctrl.ard.send_arduino_cmd_motor(0, 0)
    # if args.wait_key:
    #     try:
    #         input("Prêt à suivre. Appuyez Enter pour démarrer...")
    #     except EOFError:
    #         pass
    # if args.start_after > 0.0:
    #     time.sleep(args.start_after)

    #t0 = time.time()
    t0 = wait_for(args.start_time)

    def _on_sigint(sig, frame):
        # Sauvegarde immédiate de la trace puis arrêt moteurs
        try:
            with open(args.gpx, 'w') as f:
                f.write(gpx.to_xml())
                f.flush()
                os.fsync(f.fileno())
            print("\nCtrl+C: GPX sauvegardé -> %s" % args.gpx)
        except Exception as e:
            print("\nCtrl+C: échec sauvegarde GPX:", e)
        finally:
            try:
                ctrl.ard.send_arduino_cmd_motor(0, 0)
            except Exception:
                pass
            sys.exit(130)  # 128 + SIGINT
    signal.signal(signal.SIGINT, _on_sigint)

    print("Tracking a(t) pendant %.0f s" % args.duration)
    track_at(ctrl, gps, traj, t_start=t0, u_d_nom=args.u_d, k_pos=args.k_pos, dt_cmd=args.dt_cmd, duration_s=args.duration, seg_robot=seg_robot, seg_ref=None)

    ctrl.ard.send_arduino_cmd_motor(0, 0)

    # Trace GPX
    with open(args.gpx, 'w') as f:
        f.write(gpx.to_xml())
    print("Trace GPX -> %s" % args.gpx)

if __name__ == "__main__":
    try:
        main()
    finally:
        try:
            arddrv.ArduinoIO().send_arduino_cmd_motor(0, 0)
        except Exception:
            pass
