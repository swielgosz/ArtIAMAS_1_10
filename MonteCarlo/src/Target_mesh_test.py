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
from helper_calculations.sensor_vision import sensor_vision, sensor_reading, target_localization_both, target_localization_bearing, target_localization_radius
from helper_calculations.fisher_information import build_FIM, build_map_FIMS, plot_uncertainty_ellipse
from helper_calculations.localization_calculations import sensor_localization_routine
from helper_calculations.sensor_vision import get_LOS_coeff
from helper_calculations.penalty_fcn import min_distance_penalty, min_sensor_distance_penalty, valid_sensor_penalty
from helper_calculations.utils import get_csv_size
from helper_calculations.target_meshing import make_target_mesh
import random 

# Load data
my_path = Path(__file__).parent / "../data/terrain/GIS_terrain_100_100.csv"
rows, cols = get_csv_size(my_path)
terrain = Configuration(cols, rows) # cols = width, height = rows
terrain.load_from_csv(my_path)

# Define rectangular areas to place targets in
num_areas = 2
area_origin = [[40,40],[30,50]] # can add additional origins if using multiple rectangles
area_dim = [[10,20],[5,7]] # same as above
target_mesh_points = []

# Place targets
for i in range(num_areas):
    target_mesh_points = make_target_mesh(target_mesh_points,terrain,area_origin[i-1], area_dim[i-1])

# Print locations
print("Target locations:", target_mesh_points)

# Plot
terrain.plot_grid_with_targets(area_origin,area_dim)
