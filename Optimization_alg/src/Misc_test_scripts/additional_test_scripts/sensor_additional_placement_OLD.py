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
1. Optimally placed an additional sensor per [1] for a 
    a. homogenous sensor network localizing
    b. one target, p

TO DO (in order of importance):
0. Add functionality for multiple targets that may or may not be seen
1. Add distance sensor calculations
'''

# --------------PARAMETERS------------------
# TERRAIN INPUTS
terrain_height = 100 #796
terrain_width = 100 #1002
terrain = Configuration(terrain_width, terrain_height)
my_path = Path(__file__).parent / "../data/terrain/GIS_terrain_resize.csv"
terrain.load_from_csv(my_path)

# SENSOR LIST
sensor_rad = [25, 25, 25, 20]
sensor_type = ["acoustic","acoustic","seismic","acoustic"]
num_sensors = len(sensor_type)
sensor_comm_ratio = 1.75 # ratio of sensor communication to sensing radius 
meas_type = "bearing"
sensor_locs = [57, 58, 36, 46]
target = [50, 50]
# ------------------------------------------

# Need to first localize the target
# Calculate the bearing measurements
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

target_loc = target_localization(bearing_measurements, sensor_localization)
print("target localization:",target_loc)
print("Actual target position:", target)

# Now let's optimally add a sensor
# Only built from sensors that saw the target
FIM = build_FIM(sensor_localization, target_loc)
print("FIM:", FIM)
print("Determinant of FIM:", np.linalg.det(FIM))

# define an objective fun that optimizes a point to the sensor list
    # LATER: we'll feed this into a real optimizer
# Let's place a sensor on the map at ALL locations

def objective_fcn(x, args):
    # accept: sensor list, target
    # we're placing the LAST sensor in the list
    sensor_rad, sensor_type, num_sensors, sensor_comm_ratio, map_list, target = args

    xs, ys = x[-2], x[-1]

    valid_placement_check, valid_WSN = True, True
    
    # Newly placed sensor should not be too close to target
    if ((xs - target[0])**2 + (ys - target[1])**2)**(1/2) < 10:
        valid_placement_check = False

    # Sensor should be in the WSN
    _map = make_basic_seismic_map(num_sensors, sensor_rad, sensor_type, sensor_comm_ratio, map_list)
    valid_WSN = terrain.is_configuration_valid(_map)        

    if valid_placement_check and valid_WSN:
         FIM = build_FIM(sensor_locs, target_loc)
         return np.linalg.det(FIM)
    else:
         return 0

# Maximization!
# FOR NOW, we're just running one sensor so we can calculate all potential locations 
# and take the max. Obviously later, we'll do a formal optimization

det = [[0 for i in range(terrain_height)] for j in range(terrain_width)]
x_star = [0,0]
opt_val = 0
sensor_locs.append(1)
sensor_locs.append(1)
n = len(sensor_locs)

for i in range(terrain_height): #height == y?
    for j in range(terrain_width): #width == x?
            sensor_locs[n-2], sensor_locs[n-1] = j, i #reset the positions
            x = [j,i]
            sensor_lists = (sensor_rad, sensor_type, num_sensors, sensor_comm_ratio, sensor_locs, target) 
            obj_eval = objective_fcn(x, sensor_lists)
            print(obj_eval)
            if obj_eval > opt_val:
                x_star = x
                opt_val = obj_eval

print("Optimal det(FIM) value:", opt_val)
print("Optimal placement: ", x_star)

sensor_locs[n-2], sensor_locs[n-1] = int(x_star[0]), int(x_star[1])
print(sensor_locs)

# Make and plot a map
_map = make_basic_seismic_map(num_sensors, sensor_rad, sensor_type, sensor_comm_ratio, sensor_locs)
ax = terrain.plot_grid(_map)
new_map = add_targets(ax, target) # add the target

# Add the sensor vision stuff
for i in range(len(sensor_locs)//2):
    sensor = [0, 0]
    sensor[0], sensor[1] = sensor_locs[0+2*i], sensor_locs[1+2*i]
    nodes, _ = sensor_vision(sensor, sensor_rad[i], target, meas_type)
    
    # plot bearing or radius
    for node in nodes:
        ax.plot(node[0], node[1], "b.")

plt.title("Optimal placement of an additional sensor")
plt.show()
