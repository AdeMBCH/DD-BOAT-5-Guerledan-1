import sys
import os
import time
from pyproj import Proj, transform
import numpy as np
# access to the drivers
sys.path.append(os.path.join(os.path.dirname(__file__), 'drivers-ddboat-v2'))
import gps_driver_v2 as gpddrv
gps = gpddrv.GpsIO(tty_dev=0) # create a GPS object
gps.set_filter_speed("0") # allowing GPS measures to change even if the DDBoat is not moving
# define the projection from WGS84 (lon,lat in degrees) to UTM (meters)
# UTM zone for Brittany is 30N
projDegree2Meter = Proj("+proj=utm +zone=30, +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
# reference point "Guerledan's lake"
reference_lat = 48.199111
reference_lon = -3.014930
reference_x,reference_y = projDegree2Meter(reference_lon,reference_lat)
print ("ref",reference_x,reference_y,reference_lon,reference_lat)

import gpxpy.gpx
# Create a GPX object
gpx = gpxpy.gpx.GPX()
gpx_track = gpxpy.gpx.GPXTrack()
gpx.tracks.append(gpx_track)
gpx_segment = gpxpy.gpx.GPXTrackSegment()
gpx_track.segments.append(gpx_segment)


# convert lat,lon from DDMM.MMMM to DD.DDDDD
def cvt_gll_ddmm_2_dd (st):
    ilat = st[0]
    ilon = st[2]
    olat = float(int(ilat/100))
    olon = float(int(ilon/100))
    olat_mm = (ilat%100)/60
    olon_mm = (ilon%100)/60
    olat += olat_mm
    olon += olon_mm
    if st[3] == "W":
        olon = -olon
    return olat,olon

cnt = 200 # takes 5 GPS measures
while True:
    gll_ok,gll_data=gps.read_gll_non_blocking()
    if gll_ok: # GPGLL message received
        lat,lon = cvt_gll_ddmm_2_dd(gll_data) # convert DDMM.MMMM to DD.DDDDD
        x,y = projDegree2Meter(lon, lat) # convert to meters
        lat_check,lon_check = projDegree2Meter(x,y,inverse=True) # check conversion OK
        dx = x - reference_x
        dy = y - reference_y
        distance = np.sqrt(dx*dx+dy*dy)
        heading_trigo = np.degrees(np.arctan2(dy,dx))
        heading_geo = 90.0 - heading_trigo # convert from trigonomety to geographic
        print ("lat=%.4f lon=%.4f (check %.4f %.4f) x=%.2f y=%.2f dx=%.2f, dy=%.2f, distance=%.2f, heading=%.2f"%
        (lat,lon,lat_check,lon_check,x,y,dx,dy,distance,heading_geo))
        gpx_segment.points.append(gpxpy.gpx.GPXTrackPoint(lat, lon))
        cnt -= 1
        if cnt==0:
            with open("gps_data.gpx", "w") as f:
                f.write(gpx.to_xml())
            break
    time.sleep(0.5)
