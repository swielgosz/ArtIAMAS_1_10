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
from sensor_vision import sensor_target_angle
from fisher_information import build_FIM, build_map_FIMS, plot_uncertainty_ellipse
from localization_calculations import sensor_localization_routine


'''
This script:
1. Visualizes the "sensor vision," or what each sensor sees
2. Localizes p targets with N arbitrary sensors
3. Returns the estimates of the target locations 
4. Calculates the FIM about each point and the FIM score amongst the entire map

TO DO (in order of importance):
'''

# --------------PARAMETERS------------------
# TERRAIN INPUTS
terrain_height = 100 #796
terrain_width = 100 #1002
terrain = Configuration(terrain_width, terrain_height)
my_path = Path(__file__).parent / "../data/terrain/GIS_terrain_resize.csv"
terrain.load_from_csv(my_path)

# SENSOR LIST
sensor_rad = [15, 15, 15]
sensor_type = ["acoustic","seismic", "seismic"]
num_sensors = len(sensor_type)
sensor_comm_ratio = 1 # ratio of sensor communication to sensing radius 
meas_type = ["bearing", "radial", "radial"] #"bearing" or "radial" options
#sensor_locs = [40, 51, 54, 40, 55, 55]
targets = [46, 44]
sensor_locs = [40, 51]
phi_r1 = -44.2+49.42
phi_r2 = 44.2+49.42
sensor_rad1 = [5*np.cos(phi_r1*3.1415/180)+targets[0], 5*np.sin(phi_r1*3.1415/180)+targets[1]]
sensor_rad2 = [5*np.cos(phi_r2*3.1415/180)+targets[0], 5*np.sin(phi_r2*3.1415/180)+targets[1]]
sensor_locs.append(sensor_rad1[0])
sensor_locs.append(sensor_rad1[1])
sensor_locs.append(sensor_rad2[0])
sensor_locs.append(sensor_rad2[1])
# ------------------------------------------

# BUILD THE MAP
_map = make_basic_seismic_map(num_sensors, sensor_rad, sensor_type, meas_type, sensor_comm_ratio, sensor_locs)
score = terrain.get_configuration_observability(_map)
ax = terrain.plot_grid(_map)
new_map = add_targets(ax, targets)

## DETERMINE THE LOCALIZABLE TARGETS
# Conditions: (1) a target is seen by two bearing sensors
            # (2) a target is seen by a bearing and a radius
            # (3) a target is seen by three radii sensors

# Will return a list of lists, where each list is sensor positions
# and sensor types that can localize the target. The number of elements 
# in the larger list will correspond to the number of targets

# Coded up in "sensor_localization_routine" fcn
inputs = sensor_rad, meas_type, sensor_locs, targets, ax
(target_localized_successfully, ax) = sensor_localization_routine(inputs)

# Some information on map
# angles
theta_list = []
theta_list_offset = []
for i in range(len(sensor_locs)//2):
    sensor_in = [sensor_locs[0+2*i], sensor_locs[1+2*i]]
    target_in = targets
    theta_list.append(sensor_target_angle(sensor_in, target_in)*180/3.14)
    theta_list_offset.append(theta_list[i] - theta_list[0])

print(theta_list)
print(theta_list_offset)

## FIM CALCULATION
# For each target that is localized, calculate the FIM about that point
# For each target, create sublists of sensors that see the target

# Coded up in "fisher_information.py"
inputs = target_localized_successfully, targets, sensor_locs, sensor_rad, meas_type
(FIMs, det_sum) = build_map_FIMS(inputs)

print("MAP FIM SCORE")
print(det_sum)
print(FIMs)

new_map = plot_uncertainty_ellipse(new_map, FIMs[0], targets[0:2], 3)

plt.show()
