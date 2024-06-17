from landscape import Configuration, Node
from sensors import Sensor
import numpy as np
import pprint as pp
import os
import shutil
from pathlib import Path
from scipy.optimize import minimize, root
from scipy import optimize
import time
from make_map import make_basic_seismic_map
import matplotlib.pyplot as plt
from make_target_map import add_targets
from sensor_vision import sensor_target_angle, get_LOS_coeff
from fisher_information import build_FIM, build_map_FIMS, plot_uncertainty_ellipse
from localization_calculations import sensor_localization_routine
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning) 


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
targets = [49, 42]
sensor_locs = [52, 35, 52, 51, 41, 45]


# ------------------------------------------


inputs = sensor_rad, meas_type, sensor_locs, targets, None
(target_localized_successfully, _) = sensor_localization_routine(inputs)



## FIM CALCULATION
# Plot w/o sigma scaling
inputs = target_localized_successfully, targets, sensor_locs, sensor_rad, meas_type, terrain
(FIMs_no_scale, det_sum) = build_map_FIMS(inputs)

sensor_locs = [52, 35, 52, 51, 56, 39]
inputs = target_localized_successfully, targets, sensor_locs, sensor_rad, meas_type, terrain
(FIMs_scale, det_sum) = build_map_FIMS(inputs)

_map = make_basic_seismic_map(num_sensors, sensor_rad, sensor_type, meas_type, sensor_comm_ratio, sensor_locs)
score = terrain.get_configuration_observability(_map)
ax = terrain.plot_grid(_map)
new_map = add_targets(ax, targets)

print("MAP FIM SCORE")
print("Locations:", sensor_locs)
print("FIM det sum:", det_sum)
print("FIM itself:",FIMs_no_scale)
new_map = plot_uncertainty_ellipse(new_map, FIMs_no_scale[0], targets[0:2], 3)

# Plot w/ sigma scaling
print("MAP FIM SCORE W LOS")
print("Locations:", sensor_locs)
print("FIM det sum:", det_sum)
print("FIM itself:",FIMs_scale)
new_map = plot_uncertainty_ellipse(new_map, FIMs_scale[0], targets[0:2], 3)


plt.show()
