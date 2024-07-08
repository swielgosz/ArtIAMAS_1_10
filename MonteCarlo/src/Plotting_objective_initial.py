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
from helper_calculations.sensor_vision import sensor_vision, sensor_reading, target_localization_both, target_localization_bearing, target_localization_radius
from helper_calculations.fisher_information import build_FIM, build_map_FIMS, plot_uncertainty_ellipse
from helper_calculations.localization_calculations import sensor_localization_routine
from helper_calculations.sensor_vision import get_LOS_coeff
from helper_calculations.penalty_fcn import min_distance_penalty, min_sensor_distance_penalty, valid_sensor_penalty

# -------- DESCRIPTION -----------------
# This script plots the objective function on a contour
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
sensor_type = ["seismic","acoustic","seismic", "acoustic"]
num_sensors = len(sensor_type)
sensor_comm_ratio = 0.4 # ratio of sensor communication to sensing radius 
meas_type = ["radius", "bearing", "radius", "bearing"] #radius or bearing
LOS_flag = 0 # 1 if want to consider LOS, 0 if don't want to
sensor_locs = [40.69553033, 50.49727325, 50.22817631, 40.1569826,  44.11675794, 55.97531727, 49.70569476, 59.18030557]

# TARGETS
# Define the target mesh here!
target_mesh_pts = [[50, 50],[49, 49],[51, 51], [60, 60], [62, 62], [58, 58], [40, 40], [42, 42], [53, 48], [55, 50], [53, 50], [55, 48]]

# take out of mesh
targets = []
target_localized_successfully = []
for target in target_mesh_pts:
    targets.append(target[0])
    targets.append(target[1])
    target_localized_successfully.append(1)

# SCRIPT PARAMS
threshold = 6 #minimum distance a sensor must be from a target
d_sens_min = 3 #minimum sensor-sensor distance

# ------------------------------------------

# Generate the contour plot
print("begin contour")
size = (terrain_height, terrain_width)
print(targets)
#X = np.linspace(0, terrain_height, terrain_height)
#Y = np.linspace(0, terrain_height, terrain_height)
objective = np.zeros(size)
for i in range(terrain_height):
    for j in range(terrain_width):
        target_input = [i, j]
        inputs = target_localized_successfully, target_input, sensor_locs, sensor_rad, meas_type, terrain, 0 #LOS_flag == 1
        (FIM_target, det_sum) = build_map_FIMS(inputs)
        # j, i flip compensates for backwardness of our set-up
        objective[j, i] = np.trace(FIM_target[0]) #np.linalg.det(FIM_target[0])


fig, ax = plt.subplots()
levels = [1, 1.2, 1.4, 1.6, 1.8]
contour_plot = ax.contour(objective, levels = levels)

# Now plot the targets
tar_x = []
tar_y = []
for i in range(len(targets)//2):
    tx, ty = targets[0+2*i], targets[1+2*i]
    tar_x.append(tx)
    tar_y.append(ty)

print(tar_x, tar_y)
ax.scatter(tar_x, tar_y, marker = "^", color = "r")
print

# Finally plot the sensor positions
for i in range(len(sensor_locs)//2):
    sx, sy = sensor_locs[0+2*i], sensor_locs[1+2*i]
    if meas_type[i] in ("radial", "radius"):
        color_pick = '#000000'
    else:
        color_pick = '#e88e10'
    print(meas_type[i], color_pick)
    ax.scatter(sx, sy, color = color_pick)
    plt.text(sensor_locs[0+2*i]+1, sensor_locs[1+2*i]+1, "In")

ax.clabel(contour_plot)
plt.gca().invert_yaxis()

plt.show()
