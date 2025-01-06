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
from helper_calculations.fisher_information import build_FIM, build_map_FIMS, plot_uncertainty_ellipse, uncertainty_ellipse_stats
from helper_calculations.localization_calculations import sensor_localization_routine
from helper_calculations.target_meshing import make_target_mesh
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

# CASE 1 TARGET LIST
area_origin =  [[39, 41], [41, 45],[43, 49],[49, 51],[53, 53], [57, 57]]
area_dim = [[14, 4], [12, 4], [14, 2], [8,2], [6, 4], [6, 4]]

# CASE 2 TARGET LIST
#area_origin =  [[55, 35], [57, 33], [61, 33], [63, 33], [34, 28], [36, 32], [40, 34]]
#area_dim = [[2, 8], [4, 10], [2, 6], [2, 4], [8, 4], [8, 2], [4, 2]]


# Place targets
step = 2
target_mesh_points = []
for i in range(len(area_origin)):
    target_mesh_points = make_target_mesh(target_mesh_points,terrain,area_origin[i], area_dim[i], step)

# Break out the mesh pts:
targets = []
for target in target_mesh_points:
    for pt in target:
        targets.append(pt) 

# SENSOR LIST - CASE 1 EXPLORE
#sensor_locs = [63.62323363, 49.52390133, 42.59499122, 56.18002787, 47.2927944,  59.42022866, 58.17459166, 42.52344475,] #trace only
sensor_locs = [63.21988299009363, 50.301701179412916, 41.62841105952802, 55.718909428926395, 57.80425238524424, 42.869635611693795, 46.9996883909076, 58.33687905728023] #trace and det
sensor_rad = [50, 50, 50, 50]
sensor_type = ["seismic","seismic","seismic","acoustic"]
num_sensors = len(sensor_type)
sensor_comm_ratio = 0.4 # ratio of sensor communication to sensing radius 
meas_type = ["radius","radius", "radius","bearing"]
target_localized_successfully = [1 for _ in target_mesh_points]

# SENSOR LIST - CASE 1 EXPLOIT
#sensor_locs = [63.21988299009363, 50.301701179412916, 41.62841105952802, 55.718909428926395, 57.80425238524424, 42.869635611693795, 46.9996883909076, 58.33687905728023, 40.33305390438449, 49.568180663516735, 49.5, 39.89546817050247] #trace and det
#sensor_rad = [50, 50, 50, 50, 50, 50]
#sensor_type = ["seismic","seismic","seismic","acoustic", "acoustic", "seismic"]
#num_sensors = len(sensor_type)
#sensor_comm_ratio = 0.4 # ratio of sensor communication to sensing radius 
#meas_type = ["radius", "radius", "radius", "bearing", "bearing", "radius"]
#target_localized_successfully = [1 for _ in target_mesh_points]
#targets = [42, 41, 44, 44, 48, 48, 41, 43, 52, 46]

# SENSOR LIST - CASE 2 EXPLORE
#sensor_locs = [54.85023264, 46.29158785, 41.45345887, 39.23412295, 47.45079212, 33.72859947] #trace only
#sensor_locs = [54.50633808, 46.32458927, 40.44200985, 39.33369424, 47.27705969, 30.66311726] #trace and det
#sensor_rad = [35, 35, 35]
#sensor_type = ["seismic","acoustic","seismic"]
#num_sensors = len(sensor_type)
#sensor_comm_ratio = 0.75 # ratio of sensor communication to sensing radius 
#meas_type = ["radius", "bearing", "radius"]
#target_localized_successfully = [1 for _ in target_mesh_points]

# SENSOR LIST - CASE 2 EXPLOIT
#sensor_rad = [35, 35, 35, 35, 35]
#sensor_type = ["seismic","acoustic","seismic","seismic", "seismic"]
#num_sensors = len(sensor_type)
#sensor_comm_ratio = 0.75 # ratio of sensor communication to sensing radius 
#meas_type = ["radius", "bearing", "radius", "radius", "radius"] #radius or bearing
#sensor_locs = [40, 51, 54, 40, 55, 55]
#targets = [39, 32, 55, 40, 58, 39, 64, 33, 41, 30] 
#sensor_locs = [54.50633808, 46.32458927, 40.44200985, 39.33369424, 47.27705969, 30.66311726, 56.79309556470051, 30.492684042066756, 49.826931870141756, 40.33081847279378]
#target_localized_successfully = [1,1,1,1,1]
#LOS = 0

# ------------------------------------------

# Localize the target (not sure if calcs are needed, but sanity check)
# Coded up in "sensor_localization_routine" fcn
#inputs = sensor_rad, meas_type, sensor_locs, targets, None
#(target_localized_successfully, _) = sensor_localization_routine(inputs)


## FIM CALCULATION
inputs = target_localized_successfully, targets, sensor_locs, sensor_rad, meas_type, terrain, 0

(FIMs_exact, det_sum) = build_map_FIMS(inputs)
print("MAP FIM SCORE")
print("FIM det sum:", det_sum)
print("FIM itself:",FIMs_exact)

_map = make_basic_seismic_map(num_sensors, sensor_rad, sensor_type, meas_type, sensor_comm_ratio, sensor_locs)
score = terrain.get_configuration_observability(_map)
ax = terrain.plot_grid(_map)
new_map = add_targets(ax, targets)
tr_sum, det_sum = 0., 0.

for i, FIM_in in enumerate(FIMs_exact):
    new_map = plot_uncertainty_ellipse(new_map, FIM_in, [targets[0+2*i], targets[1+2*i]], 1, 1, "black", "analytical")
    tr_sum += np.trace(FIM_in)
    det_sum += np.linalg.det(FIM_in)

print('Trace and det sums:', tr_sum, det_sum)
colors = ["grey"]
FIM_areas, maj_axis, std_devs = [], [], []
FIM_save_data1, FIM_save_data2, target_save_data = [], [], []

for i in range(len(targets)//2):
        target_i = [targets[0+2*i], targets[1+2*i]]
        new_map = plot_uncertainty_ellipse(new_map, FIMs_exact[i], target_i, 1, 1, colors[0], "analytical")
        #new_map = plot_uncertainty_ellipse(new_map, FIMs_random[i], target_i, 2.48, 1, "blue", "numerical")

        # Compute the stats
        Cov_Mat, a, b, sx, sy, area_first = uncertainty_ellipse_stats(FIMs_exact[i], target, 2.48, 1)
        FIM_areas.append(area_first)
        maj_axis.append(a)
        std_devs.append(sx**2)
        std_devs.append(sy**2)


print("mean area:", np.mean(FIM_areas))
print("mean axis", np.max(maj_axis))
print("variances both x, y", np.mean(std_devs))

plt.show()
