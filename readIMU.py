#!/usr/bin/env python
import os, sys, time

# Driver IMU local
sys.path.append(os.path.join(os.path.dirname(__file__), 'drivers-ddboat-v2'))
import imu9_driver_v2 as imudrv

DURATION_SEC = 20.0
SLEEP_SEC    = 0.01/2  # ~100 Hz

def main():
    imu = imudrv.Imu9IO()
    t_end = time.time() + DURATION_SEC
    n = 0
    f = open("datacalibration.txt", "w")
    try:
        f.write("xmag,ymag,zmag,xaccel,yaccel,zaccel,xgyro,ygyro,zgyro\n")
        while time.time() < t_end:
            mx, my, mz = imu.read_mag_raw()
            ax, ay, az = imu.read_accel_raw()
            gx, gy, gz = imu.read_gyro_raw()
            line = "{},{},{},{},{},{},{},{},{}\n".format(mx, my, mz, ax, ay, az, gx, gy, gz)
            f.write(line)
            n += 1
            time.sleep(SLEEP_SEC)
    finally:
        f.close()
    print("echantillons ecrits: {}".format(n))

if __name__ == "__main__":
    main()
