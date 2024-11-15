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
from helper_calculations.fisher_information import build_FIM

# -------- DESCRIPTION -----------------
# This script presents the "sensor vision" of a set of sensors and
# targets on the map. It will also detail, given a set of sensors and targets, 
# which targets are localized (by which sensors).

# Note that this script isn't particularly necessary anymore - most of the
# functionality of the script was coded up in helper functions instead
# so we don't really need this anymore. 
# --------------------------------------


# --------------PARAMETERS------------------
# TERRAIN INPUTS
terrain_height = 100 #796
terrain_width = 100 #1002
terrain = Configuration(terrain_width, terrain_height)
my_path = Path(__file__).parent / "../data/terrain/GIS_terrain_resize.csv"
terrain.load_from_csv(my_path)

# SENSOR LIST
sensor_rad = [25, 25, 25, 20, 20]
sensor_type = ["acoustic","seismic","seismic","acoustic", "seismic"]
num_sensors = len(sensor_type)
sensor_comm_ratio = 1 # ratio of sensor communication to sensing radius 
meas_type = ["bearing", "radial", "radial", "bearing", "radial"] #"bearing" or "radial" options
sensor_locs = [57, 58, 36, 46, 66, 36, 75, 64, 44, 16]
targets = [40, 60, 50, 50, 63, 75, 48, 33]
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

# Check every target
target_local_sensors = []
target_local_rads = []
target_local_sensors_types = []
localization_types = []
target_localized_successfully = []
for i in range(len(targets)//2):
    # Create a list to check if the target it seen by a sensor
    # List index corresponds to sensor index: 1 = seen by that sensor, 0 = not
    bearing_sensor_list = []
    radial_sensor_list = []
    sensor_radius_list = []
    bearing_count, radial_count = 0, 0
    
    # Then check every sensor to see if the target is seen
    xt, yt = targets[0+2*i], targets[1+2*i]
    for j in range(len(sensor_locs)//2):
        x, y = sensor_locs[0+2*j], sensor_locs[1+2*j]
        dist = ((x-xt)**2+(y-yt)**2)**(1/2)
        
        if dist < sensor_rad[j]:
            if meas_type[j] == "bearing":
                bearing_count += 1
                bearing_sensor_list.append(x)
                bearing_sensor_list.append(y)
                sensor_radius_list.append(sensor_rad[j])
            elif meas_type[j] == "radial":
                radial_count += 1
                radial_sensor_list.append(x)
                radial_sensor_list.append(y)
                sensor_radius_list.append(sensor_rad[j])
    
        CND1, CND2, CND3 = bearing_count >= 2, radial_count >= 3, (bearing_count >= 1 & radial_count >= 1)
        if CND1:
            sensor_sublist = bearing_sensor_list
            sensor_type_sublist = ["bearing","bearing"]
            localization_type = "bearing"
            successful_localization = 1
            break
        if CND2:
            sensor_sublist = radial_sensor_list
            sensor_type_sublist = ["radial","radial","radial"]
            localization_type = "radial"
            successful_localization = 1
            break
        if CND3:
            sensor_sublist = bearing_sensor_list[0:2] + radial_sensor_list[0:2] 
            sensor_type_sublist = ["bearing","radial"]
            localization_type = "both"
            successful_localization = 1
            break
    
    # If, after the entire process we cannot localize a point, return empty lists
    if (CND1|CND2|CND3) == False:
        sensor_sublist = []
        sensor_type_sublist = []
        localization_type = "none"
        successful_localization = 0

    # Finally, add to the primary list of sensors that localize ith target
    target_local_sensors.append(sensor_sublist)
    target_local_sensors_types.append(sensor_type_sublist)
    localization_types.append(localization_type)
    target_local_rads.append(sensor_radius_list)
    target_localized_successfully.append(successful_localization)

print("SUBLISTS:")
print(target_local_sensors_types)
print(target_local_sensors)
print(localization_types)
print(target_localized_successfully)

## LOCALIZE EACH TARGET
# Use the sublists associated with the i'th target from above
# To localize/visualize each one with the associated calculations

print("TARGET LOCALIZATION:")
target = [0., 0.]
for i in range(len(targets)//2):
    # Define the target and it's localizing sensors
    target[0], target[1] = targets[0+2*i], targets[1+2*i]
    sub_sensor = target_local_sensors[i]
    sub_meas_type = target_local_sensors_types[i]
    sub_sensor_rads = target_local_rads[i]

    if localization_types[i] == "none":
        print("Cannot localize target",i+1,"at",target)
    else:
        print("Localizing target",i+1,"using",localization_types[i])
        measurements = []
        sensor_localization = []
        for j in range(len(sub_sensor)//2):
            
            # determine what nodes will be seen from a sensor when a target is detected
            sensor = [0, 0]
            sensor[0], sensor[1] = sub_sensor[0+2*j], sub_sensor[1+2*j]
            nodes, _ = sensor_vision(sensor, sub_sensor_rads[j], target, sub_meas_type[j])

            # Generate a list of unit vectors
            # These are passed into the localization calculation later
            measurement, seen = sensor_reading(sensor, sub_sensor_rads[j], target, sub_meas_type[j])
            if seen:
                measurements.append(measurement)
            
            # plot bearing or radius
            for node in nodes:
                ax.plot(node[0], node[1], "b.")

        # Use a different localization routine depending on type 
        if localization_types[i] == "bearing":
            target_loc = target_localization_bearing(measurements, sub_sensor)
        elif localization_types[i] == "radial":
            target_loc = target_localization_radius(measurements, sub_sensor)
        elif localization_types[i] == "both":
            target_loc = target_localization_both(measurements, sub_sensor, sub_meas_type)
        
        print("target localization:",target_loc)
        print("Actual target position:", target)

plt.show()
