from helper_calculations.landscape import Configuration, Node
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
# This script was used to compute the exact solution for the simple
# case where there is one target, two range, and one bearing sensors
# This was then fed into the simple_figure_plotting_for_paper script
# --------------------------------------

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
targets = [49, 42]
sensor_locs = [52, 35]
LOS = 0

# ------- CALC EXACT SOLT --------
rb = ((sensor_locs[0] - targets[0])**2 + (sensor_locs[1] - targets[1])**2)**(1/2)
phi_b = sensor_target_angle(sensor_locs[0:2], targets)
alpha1 = 1/rb**2
alpha2 = 1

phi_r1 = (-1/2)*np.arccos(alpha1/(2*alpha2))
if np.abs(phi_r1) < 1.5:
    phi_r1-=np.pi

print(phi_r1)
#sensor_locs.append(rb*np.cos(phi_r1+phi_b)+targets[0])
#sensor_locs.append(rb*np.sin(phi_r1+phi_b)+targets[1])
sensor_locs.append(52)
sensor_locs.append(51)

# ------------------------------------------

# Localize the target (not sure if calcs are needed, but sanity check)
# Coded up in "sensor_localization_routine" fcn
inputs = sensor_rad, meas_type, sensor_locs, targets, None
(target_localized_successfully, _) = sensor_localization_routine(inputs)

def fun(x):
    return -(1/rb**2)*np.sin(2*(x)) + np.sin(2*(x-phi_r1))

def second_deriv(x):
    return (1/rb**2)*np.cos(2*(phi_b-x))*(-2) + np.cos(2*(x-phi_r1))*2

def full_det(x):
    return (1/rb**2)*np.cos((phi_b-phi_r1)) + (1/rb**2)*np.cos((phi_b-x))+ np.sin(2*(x-phi_r1))

sol = optimize.bisect(fun, 0, 3.14/2)
sol = sol
if np.abs(phi_r1) > 1.5:
    sol+=np.pi
print("optimized angle:", sol)
print("bearing radius dist, rb", rb)
print("Eval of det w/ angle:", full_det(sol))
print("Eval of second deriv w/ angle:", second_deriv(sol))

# Now place this sensor on the map 
sensor_locs.append(float((rb+1)*np.cos(sol+phi_b)+targets[0]))
sensor_locs.append(float((rb+1)*np.sin(sol+phi_b)+targets[1]))

## FIM CALCULATION
inputs = target_localized_successfully, targets, sensor_locs, sensor_rad, meas_type, terrain, LOS
(FIMs_exact, det_sum) = build_map_FIMS(inputs)

print("MAP FIM SCORE")
print("Locations:", sensor_locs)
print("FIM det sum:", det_sum)
print("FIM itself:",FIMs_exact)

# FIM w/ DIRECT optimal value
#sensor_locs = [52, 35, 52, 51, 40.129629629629626, 45.833333333333336]
#inputs = target_localized_successfully, targets, sensor_locs, sensor_rad, meas_type, terrain, LOS
#(FIMs_numerical, det_sum) = build_map_FIMS(inputs)
#print("MAP FIM SCORE")
#print("Locations:", sensor_locs)
#print("FIM det sum DIRECT:", det_sum)
#print("FIM itself:",FIMs_numerical)

# Redo, but with some other angles and see what happens
# Perturn by 30 deg
#sensor_locs[4] = float((rb+1)*np.cos(sol-0+phi_b)+targets[0])
#sensor_locs[5] = float((rb+1)*np.sin(sol-0+phi_b)+targets[1])

# (40.6, 36.7)
sensor_locs[4] = 40.6
sensor_locs[5] = 36.7
inputs = target_localized_successfully, targets, sensor_locs, sensor_rad, meas_type, terrain, LOS
(FIMs_1, det_sum) = build_map_FIMS(inputs)
print("FIM 1:", FIMs_1, det_sum)

# (60.5, 43)
sensor_locs[4] = 60.5
sensor_locs[5] = 43
inputs = target_localized_successfully, targets, sensor_locs, sensor_rad, meas_type, terrain, LOS
(FIMs_2, det_sum) = build_map_FIMS(inputs)
print("FIM 2:", FIMs_2, det_sum)

# (47.1, 52.6)
sensor_locs[4] = 47.1
sensor_locs[5] = 52.6
inputs = target_localized_successfully, targets, sensor_locs, sensor_rad, meas_type, terrain, LOS
(FIMs_3, det_sum) = build_map_FIMS(inputs)
print("FIM 3:", FIMs_3, det_sum)


_map = make_basic_seismic_map(num_sensors, sensor_rad, sensor_type, meas_type, sensor_comm_ratio, sensor_locs)
score = terrain.get_configuration_observability(_map)
ax = terrain.plot_grid(_map)
new_map = add_targets(ax, targets)

new_map = plot_uncertainty_ellipse(new_map, FIMs_exact[0], targets[0:2], 2.89, 1, "black", "analytical")
new_map = plot_uncertainty_ellipse(new_map, FIMs_1[0], targets[0:2], 2.89, 1, "tab:blue", "analytical")
new_map = plot_uncertainty_ellipse(new_map, FIMs_2[0], targets[0:2], 2.89, 1, "tab:red", "analytical")
new_map = plot_uncertainty_ellipse(new_map, FIMs_3[0], targets[0:2], 2.89, 1, "tab:purple", "analytical")


plt.show()
