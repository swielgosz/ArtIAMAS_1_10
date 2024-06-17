from helper_calculations.landscape import Configuration
from helper_calculations.sensors import Sensor
import numpy as np
import pprint as pp
import os
import shutil
from pathlib import Path
from scipy.optimize import minimize, root
from scipy import optimize
import time
from helper_calculations.make_map import make_basic_seismic_map
import matplotlib.pyplot as plt
from helper_calculations.make_target_map import add_targets
from helper_calculations.sensor_vision import sensor_target_angle
from helper_calculations.fisher_information import build_FIM, build_map_FIMS, plot_uncertainty_ellipse
from helper_calculations.localization_calculations import sensor_localization_routine
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning) 

# -------- DESCRIPTION -----------------
# This script was makes plots the sensors (radii included) and 
# the uncertainty ellipses. It's mostly just for inspection purposes
# --------------------------------------

# --------------PARAMETERS------------------
# TERRAIN INPUTS
terrain_height = 100 #796
terrain_width = 100 #1002
terrain = Configuration(terrain_width, terrain_height)
my_path = Path(__file__).parent / "../data/terrain/GIS_terrain_resize.csv"
terrain.load_from_csv(my_path)

# SENSOR LIST
sensor_rad = [25, 25, 25, 15, 15, 20]
sensor_type = ["seismic","acoustic","seismic","seismic","acoustic","seismic"]
num_sensors = len(sensor_type)
sensor_comm_ratio = 2 # ratio of sensor communication to sensing radius 
meas_type = ["radius", "bearing", "radius","radius", "bearing", "radius"] #radius or bearing
#sensor_locs = [40, 51, 54, 40, 55, 55]
targets = [37, 39, 42, 32, 47, 39, 41, 45] 
sensor_locs = [29, 44, 33, 29, 58, 34, 27.5, 40.9, 44.6, 58.1, 48.2, 23.8]
LOS = 1

# ------------------------------------------

# Localize the target (not sure if calcs are needed, but sanity check)
# Coded up in "sensor_localization_routine" fcn
inputs = sensor_rad, meas_type, sensor_locs, targets, None
(target_localized_successfully, _) = sensor_localization_routine(inputs)
target_localized_successfully = [1,1,1,1]

## FIM CALCULATION
inputs = target_localized_successfully, targets, sensor_locs, sensor_rad, meas_type, terrain, LOS

(FIMs_exact, det_sum) = build_map_FIMS(inputs)
print("MAP FIM SCORE")
print("FIM det sum:", det_sum)
print("FIM itself:",FIMs_exact)

_map = make_basic_seismic_map(num_sensors, sensor_rad, sensor_type, meas_type, sensor_comm_ratio, sensor_locs)
score = terrain.get_configuration_observability(_map)
ax = terrain.plot_grid(_map)
new_map = add_targets(ax, targets)

for i, FIM_in in enumerate(FIMs_exact):
    new_map = plot_uncertainty_ellipse(new_map, FIM_in, [targets[0+2*i], targets[1+2*i]], 2.48, 1, "black", "analytical")



plt.show()
