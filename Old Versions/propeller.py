import sys
import os
import time
# access to the drivers
sys.path.append(os.path.join(os.path.dirname(__file__), 'drivers-ddboat-v2'))
import arduino_driver_v2 as arddrv
ard = arddrv.ArduinoIO() # create an ARduino object
left_speed = 100
right_speed = 50#-left_speed
ard.send_arduino_cmd_motor(left_speed,right_speed) # in place turn
time.sleep(2)
ard.send_arduino_cmd_motor(0,0) # stop