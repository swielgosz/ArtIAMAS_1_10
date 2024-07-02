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
sensor_rad = [100, 100, 100]
sensor_type = ["seismic","acoustic","seismic"]
num_sensors = len(sensor_type)
sensor_comm_ratio = 1.5 # ratio of sensor communication to sensing radius 
meas_type = ["radius", "bearing", "radius"] #radius or bearing
LOS_flag = 1 # 1 if want to consider LOS, 0 if don't want to

# Individal sensors
targets = []
sensor_locs = []
tar1 = [37, 39]
tar2 = [42, 32] 
tar3 = [47, 39] 
tar4 = [41, 45] 
targets.extend(tar1)
targets.extend(tar2)
targets.extend(tar3)
targets.extend(tar4)
target_localized_successfully = [1,1,1,1]

sens1 = [29, 44] 
sens2 = [33, 29] 
sens3 = [58, 34] 
#sens4 = [44, 29] 
sensor_locs.extend(sens1)
sensor_locs.extend(sens2)
sensor_locs.extend(sens3)
#sensor_locs.extend(sens4)


# SCRIPT PARAMS
threshold = 6 #minimum distance a sensor must be from a target
d_sens_min = 3 #minimum sensor-sensor distance

# ------------------------------------------

# ADDED SENSOR LISTS (from optimization)
sensor_rad_new = [100, 100, 100]
sensor_type_new = ["seismic", "seismic", "acoustic"]
num_sensors_new = len(sensor_type_new)
sensor_comm_ratio_new = sensor_comm_ratio # ratio of sensor communication to sensing radius 
meas_type_new = ["radial", "radial","bearing"]
sensors_new = [48.50243865, 50.51935681, 46.39833867, 24.44360616, 47.66750495, 38.5]
sensors_final = sensor_locs.copy()
sensor_rad_total = sensor_rad.copy()
meas_type_total = meas_type.copy()

for i in range(len(sensors_new)//2):
    sx, sy = sensors_new[0+2*i], sensors_new[1+2*i]
    sensors_final.append(sx)
    sensors_final.append(sy)
    sensor_rad_total.append(sensor_rad_new[i])
    meas_type_total.append(meas_type_new[i])

# Generate the contour plot
print("begin contour")
size = (terrain_height, terrain_width)
print(sensors_final)
print(targets)
#X = np.linspace(0, terrain_height, terrain_height)
#Y = np.linspace(0, terrain_height, terrain_height)
objective = np.zeros(size)
for i in range(terrain_height):
    for j in range(terrain_width):
        target_input = [i, j]
        inputs = target_localized_successfully, target_input, sensors_final, sensor_rad_total, meas_type_total, terrain, 0 #LOS_flag == 1
        (FIMs_no_LOS, det_sum) = build_map_FIMS(inputs)
        objective[i,j] = det_sum

fig, ax = plt.subplots()
#levels = [0, 1, 2, 3, 4, 5]
contour_plot = ax.contour(objective)#, levels = levels)

# Now plot the targets
tar_x = []
tar_y = []
for i in range(len(targets)//2):
    tx, ty = targets[0+2*i], targets[1+2*i]
    tar_x.append(tx)
    tar_y.append(ty)

print(tar_x, tar_y)
ax.scatter(tar_x, tar_y, marker = "^", color = "r")

# Finally plot the sensor positions
for i in range(len(sensor_locs)//2):
    sx, sy = sensor_locs[0+2*i], sensor_locs[1+2*i]
    if meas_type[i] in ("radial", "radius"):
        color_pick = '#000000'
    else:
        color_pick = '#e88e10'
    print(meas_type_new[i], color_pick)
    ax.scatter(sx, sy, color = color_pick)
    plt.text(sensor_locs[0+2*i]+1, sensor_locs[1+2*i]+1, "In")

for i in range(len(sensors_new)//2):
    sx, sy = sensors_new[0+2*i], sensors_new[1+2*i]
    if meas_type_new[i] in ("radial", "radius"):
        color_pick = '#000000'
    else:
        color_pick = '#e88e10'
    print(meas_type_new[i], color_pick)
    
    ax.scatter(sx, sy, color = color_pick)
    plt.text(sensors_new[0+2*i]+1, sensors_new[1+2*i]+1, "F")

ax.clabel(contour_plot)
plt.gca().invert_yaxis()

plt.show()
