from landscape import Configuration, Node
from sensors import Sensor
import numpy as np
import pprint as pp
import os
import shutil
from pathlib import Path
from scipy.optimize import minimize
from scipy import optimize
import time
from make_map import make_basic_seismic_map
import matplotlib.pyplot as plt
from make_target_map import add_targets
from sensor_vision import sensor_vision, sensor_reading, target_localization
from fisher_information import build_FIM


'''
This script:
1. Visualizes the "sensor vision," or what each sensor sees
2. Localizes targets with 2 sensors (bearing only now)
3. Returns the estimates of the target locations (bearing only now)

TO DO (in order of importance):
0. Add functionality for multiple targets that may or may not be seen
1. Add distance sensor vision
2. Localization calculations with bearing/distance sensor combos
3. Localization with 2+ distance sensors
4. Script control to tell you if a target is detected, can it then be localized
'''

# --------------PARAMETERS------------------
# TERRAIN INPUTS
terrain_height = 100 #796
terrain_width = 100 #1002
terrain = Configuration(terrain_width, terrain_height)
my_path = Path(__file__).parent / "../data/terrain/GIS_terrain_resize.csv"
terrain.load_from_csv(my_path)

# SENSOR LIST
sensor_rad = [25, 25, 10, 10]
sensor_type = ["acoustic","acoustic","seismic","acoustic"]
num_sensors = len(sensor_type)
sensor_comm_ratio = 1.75 # ratio of sensor communication to sensing radius 
meas_type = "radial" #"bearing" or "radial" options
# ------------------------------------------
       
sensor_locs = [57, 58, 36, 46, 65, 45]
target = [50, 50]

# Create a map with static targets
_map = make_basic_seismic_map(num_sensors, sensor_rad, sensor_type, sensor_comm_ratio, sensor_locs)
score = terrain.get_configuration_observability(_map)
ax = terrain.plot_grid(_map)

# Plots the "sensor vision"
bearing_measurements = []
sensor_localization = []
for i in range(len(sensor_locs)//2):
    
    # determine what nodes will be seen from a sensor when a target is detected
    sensor = [0, 0]
    sensor[0], sensor[1] = sensor_locs[0+2*i], sensor_locs[1+2*i]
    nodes, _ = sensor_vision(sensor, sensor_rad[i], target, meas_type)

    # Generate a list of unit vectors
    # These are passed into the localization calculation later
    bearing, seen = sensor_reading(sensor, sensor_rad[i], target, meas_type)
    if seen:
        bearing_measurements.append(bearing)
        sensor_localization.append(sensor[0])
        sensor_localization.append(sensor[1])
        print("bearing vectors",bearing_measurements)
    
    # plot bearing or radius
    for node in nodes:
        ax.plot(node[0], node[1], "b.")
    
# Adds list of targets passed to it
# add_target(plot figure, list of targets)
new_map = add_targets(ax, target)
ax.set_title("Sensor vision of a target - bearing only information")

# Now try localizing a point
target_loc = target_localization(bearing_measurements, sensor_localization)
print("target localization:",target_loc)
print("Actual target position:", target)

# Great! We can localize a point with 2 sensors
plt.show()
