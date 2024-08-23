from helper_calculations.landscape import Configuration, Node
from helper_calculations.sensors import Sensor
import numpy as np
import pprint as pp
import os
import shutil
from pathlib import Path
from scipy.optimize import minimize
from scipy import optimize
import time
from helper_calculations.make_map import make_basic_seismic_map
import matplotlib.pyplot as plt
from helper_calculations.make_target_map import add_targets
# from helper_calculations.sensor_vision import sensor_vision, sensor_reading, target_localization_both, target_localization_bearing, target_localization_radius
from helper_calculations.fisher_information import build_FIM, build_map_FIMS, plot_uncertainty_ellipse, return_obj_vals, uncertainty_ellipse_stats
from helper_calculations.localization_calculations import sensor_localization_routine
from helper_calculations.sensor_vision import get_LOS_coeff
from helper_calculations.penalty_fcn import min_distance_penalty, min_sensor_distance_penalty, valid_sensor_penalty, WSN_penalty
from helper_calculations.save_vals import save_data
from helper_calculations.target_meshing import make_target_mesh

# -------- DESCRIPTION -----------------
# This script was makes an INITIAL PLACEMENT
# --------------------------------------

# --------------PARAMETERS------------------
# TERRAIN INPUTS
terrain_height = 100 #796
terrain_width = 100 #1002
terrain = Configuration(terrain_width, terrain_height)
my_path = Path(__file__).parent / "../data/terrain/GIS_terrain_resize.csv"
terrain.load_from_csv(my_path)

# INITIAL SENSOR LIST
# Characteristic list
sensor_rad = [50, 50, 50, 50]
sensor_type = ["seismic","seismic","seismic", "acoustic"]
num_sensors = len(sensor_type)
sensor_comm_ratio = 0.4 # ratio of sensor communication to sensing radius 
meas_type = ["radius", "radius", "radius", "bearing"] #radius or bearing
LOS_flag = 0 # 1 if want to consider LOS, 0 if don't want to

# OPTIMIZATION PARAMETERS
threshold = 5 #minimum distance a sensor must be from a target
d_sens_min = 4 #minimum sensor-sensor distance
printer_counts = 500 #print the results every ___ fcn calls

# DIRECT Parameters
vol_tol = 1e-30
max_iters = 15000
maxfun_1 = 15000

# DE Parameters - optimized!
popsize = 7
mutation = 0.8
recombination = 0.8

# Sim Annealing Params
initial_temp = 3000
restart_temp_ratio = 0.0005
visit = 2.75

# ------------------------------------------
# Define rectangular areas to place targets in
area_origin =  [[55, 35], [57, 33], [61, 33], [63, 33], [34, 28], [36, 32], [40, 34]] # can add additional origins if using multiple rectangles
area_dim = [[2, 8], [4, 10], [2, 6], [2, 4], [8, 4], [8, 2], [4, 2]] # [width, height]
target_mesh_points = []

# Place targets
for i in range(len(area_origin)):
    target_mesh_points = make_target_mesh(target_mesh_points,terrain,area_origin[i], area_dim[i], 2)

# Localize the new map
target_localized_successfully = [1 for _ in range(len(target_mesh_points ))] #just need a list of 1's 
targets_in = []
for target in target_mesh_points:
    for pt in target:
        targets_in.append(pt) 

_map = make_basic_seismic_map(num_sensors, sensor_rad, sensor_type, meas_type, sensor_comm_ratio, [])
ax = terrain.plot_grid(_map)
new_map = add_targets(ax, targets_in)
plt.show()