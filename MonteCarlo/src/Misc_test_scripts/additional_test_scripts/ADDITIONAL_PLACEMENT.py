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
from sensor_vision import sensor_vision, sensor_reading, target_localization_both, target_localization_bearing, target_localization_radius
from fisher_information import build_FIM, build_map_FIMS
from localization_calculations import sensor_localization_routine



# --------------PARAMETERS------------------
# TERRAIN INPUTS
terrain_height = 100 #796
terrain_width = 100 #1002
terrain = Configuration(terrain_width, terrain_height)
my_path = Path(__file__).parent / "../data/terrain/GIS_terrain_resize.csv"
terrain.load_from_csv(my_path)

# INITIAL SENSOR LIST
sensor_rad = [15, 15, 15]
sensor_type = ["acoustic","seismic"]
num_sensors = len(sensor_type)
sensor_comm_ratio = 1.75 # ratio of sensor communication to sensing radius 
meas_type = ["bearing", "radial"] #radius or bearing
targets = [49, 42, 45, 45]
sensor_locs = [52, 35, 52, 51]

# OPTIMIZATION PARAMETERS
threshold = 8 #minimum distance a sensor must be from a target
vol_tol = 1e-20
maxfun = 50000
max_iters = 5000
printer_counts = 1000 #print the results every ___ fcn calls
# Use list to run a parameter sweep
vol_tol = [1e-20]
# vol_tol = [1e-14] #optimization param
# ------------------------------------------

## DETERMINE THE LOCALIZABLE TARGETS
inputs = sensor_rad, meas_type, sensor_locs, targets, None
(target_localized_successfully, _) = sensor_localization_routine(inputs)


## NOMINAL FIM MAP SCORE
inputs = target_localized_successfully, targets, sensor_locs, sensor_rad,meas_type
(FIMs, det_sum) = build_map_FIMS(inputs)
print("Nominal total FIM score")
# print(FIMs)
print(det_sum)

# ------------------------------------------
## ADDITIONAL PLACEMENT:
# Using a list of available sensors, this section will 
# Numerically solve for the optimal placement of the sensor
# by maximizing the FIM determinant

# Now define the objective fcn
def objective_fcn(x, *args):
    global counter
    global fcn_eval_list
    global fcn_counter

    counter += 1
    # (0) PRELIMINARIES
    # Initialize the sensor lists
    sensor_rad, sensor_type, meas_type, sensor_positions, num_sensors = [], [], [], [], 0

    # Pull out the old and new sensor lists
    existing_sensor_lists, new_sensor_lists, target_inputs = args

    sensor_positions_old, sensor_rad_old, sensor_type_old, num_sensors_old, sensor_comm_ratio, meas_type_old = existing_sensor_lists
    sensor_rad_new, sensor_type_new, num_sensors_new, sensor_comm_ratio_new, meas_type_new = new_sensor_lists
    targets, target_localized_successfully = target_inputs

    # Augment the lists with the placed and additionally placed sensors
    sensor_rad.extend(sensor_rad_old)
    sensor_rad.extend(sensor_rad_new)
    
    sensor_type.extend(sensor_type_old)
    sensor_type.extend(sensor_type_new)
    
    meas_type.extend(meas_type_old)
    meas_type.extend(meas_type_new)
    
    sensor_positions.extend(sensor_positions_old)
    
    num_sensors += num_sensors_old
    num_sensors += num_sensors_new

    # Add the optimization variables to the lists
    for k in range(len(x)//2):
        i, j = int(x[0+2*k]), int(x[1+2*k])
        sensor_positions.extend([i, j])

    # Create a map with the existing sensors and new sensor
    _map = make_basic_seismic_map(num_sensors, sensor_rad, sensor_type, meas_type, sensor_comm_ratio, sensor_positions)

    # (1) CHECK FOR VALID CONFIGURATIONS
    valid_placement_check, valid_WSN = True, True

    # (1.1) First check for valid placement on additionally placed sensors
    for check in range(len(x)//2):
        i, j = int(x[0+2*check]), int(x[1+2*check])

        # If ANY sensors are invalid, the entire config is invalid
        # Note that (i,j) are switched to respect dims of terrain class
        # this script is (x,y) = width, height. 
        # valid sensor check is (x, y) = height, width
        if not terrain.is_valid_sensor_location(j, i):
            valid_placement_check = False
    
    # (1.2) Now check for WSN
    valid_WSN = terrain.is_configuration_valid(_map)

    # (1.3) Ensure the bearing targets are a minimum distance away
    # Otherwise the 1/r^2 will cause it to just blow up
    # Sweep over all the additionally-placed sensors
    for check in range(len(x)//2):
        i, j = int(x[0+2*check]), int(x[1+2*check])

        # Compute the distance from it to a target
        for k in range(len(targets)//2):
            tx, ty = targets[0+2*k], targets[1+2*k]
            dist = ((tx-i)**2+(ty-j)**2)**(1/2)
            # If w/in the limit, this placement is invalid!
            if dist < threshold:
                valid_placement_check = False

    # (2) CALCULATE THE SCORE
    # (2.1) First see if we've localized any new targets
        # I take it back, we shouldn't be doing this
        # We only localize a target if we actually place down the sensor
        # Thus, we don't know that a priori when we're picking candidate placements
        # in the optimizaiton routine!

    #inputs = sensor_rad, meas_type, sensor_positions, targets, None
    #(target_localized_successfully, _) = sensor_localization_routine(inputs)

    # (2.2) Construct the FIM's per target + calculate the det score
    #sensor_positions[-2] = x[0]
    #sensor_positions[-1] = x[1]

    inputs = target_localized_successfully, targets, sensor_positions, sensor_rad, meas_type
    (FIMs, det_sum) = build_map_FIMS(inputs)
    det_mult = 1
    for kk in range(len(FIMs)):
        det_mult = det_mult*np.linalg.det(FIMs[kk])

    # If a valid config, return (-1)det(FIMs)
    if valid_placement_check and valid_WSN:
        # Maximize the determinant of the map
        if counter % printer_counts == 0:
            print(counter, det_mult, sensor_positions)

        fcn_eval_list.append(det_mult)
        fcn_counter.append(counter)
        return -det_mult # Minimize the negative det(FIMS)
    
    # If an invalid construction, then fail and return 0
    else:
        fcn_eval_list.append(0)
        fcn_counter.append(counter)
        if counter % printer_counts == 0:
            print(counter, 0, sensor_positions)
        return 0


# ------------------------------------------
# Set a list of additional sensors to place down
# SENSOR LIST
sensor_rad_new = [15]
sensor_type_new = ["seismic"]
num_sensors_new = len(sensor_type_new)
sensor_comm_ratio_new = [1.75] # ratio of sensor communication to sensing radius 
meas_type_new = ["radial"]

# Construct tuples to pass in
existing_sensor_lists = (sensor_locs, sensor_rad, sensor_type, num_sensors, sensor_comm_ratio, meas_type)
new_sensor_lists = (sensor_rad_new, sensor_type_new, num_sensors_new, sensor_comm_ratio_new, meas_type_new)
target_ins = (targets, target_localized_successfully)

# Set the bounds (controls the number of variables to reason over)
bounds = [(0, terrain_width-1),(0, terrain_height-1)]*(num_sensors_new) #bounds placing points to only on the map
sensor_list = (existing_sensor_lists, new_sensor_lists, target_ins)   

# Call the optimizer!
counter = 0

# This part of the script just plots the obj fcn vs count #
# so we can check progress
fig = plt.figure()
ax_opt = fig.add_subplot()
for count, vol_tol_test in enumerate(vol_tol):
    fcn_eval_list, fcn_counter = [], []
    print("Start optimizer")
    res = optimize.direct(objective_fcn, bounds=bounds, args = sensor_list, vol_tol=vol_tol_test, maxfun = maxfun, maxiter = max_iters)
    ax_opt.plot(fcn_counter, fcn_eval_list, linestyle='None', marker='+')
    print(res.message, "Final values:", res.x)
    x_out = res.x

# ------------------------------------------
# Formatting
ax_opt.set_ylabel('Objective Function Evaluations'), ax_opt.set_xlabel('Iteration Count')
ax_opt.set_ylim([-1, max(fcn_eval_list)*1.25])
ax_opt.set_title('Function Evaluation across Optimization')
ax_opt.grid(True, which='minor')  
plt.show()

# Create the finalized lists and create a new map
num_sens_total = num_sensors+num_sensors_new
sensor_rad_total, sensor_type_total, meas_type_total, sensor_positions_final = [], [], [], []
for i in range(num_sensors):
    sensor_rad_total.append(sensor_rad[i])
    sensor_type_total.append(sensor_type[i])
    meas_type_total.append(meas_type[i])
    sensor_positions_final.append(sensor_locs[0+2*i])
    sensor_positions_final.append(sensor_locs[1+2*i])

for i in range(num_sensors_new):
    sensor_rad_total.append(sensor_rad_new[i])
    sensor_type_total.append(sensor_type_new[i])
    meas_type_total.append(meas_type_new[i])
    sensor_positions_final.append(x_out[0+2*i])
    sensor_positions_final.append(x_out[1+2*i])


# Finally, make the map for plotting
_map = make_basic_seismic_map(num_sens_total, sensor_rad_total, sensor_type_total, meas_type_total, sensor_comm_ratio, sensor_positions_final)
ax = terrain.plot_grid(_map)
# Localize the new map
inputs = target_localized_successfully, targets, sensor_positions_final, sensor_rad_total, meas_type_total
(FIMs, det_sum) = build_map_FIMS(inputs)
print(targets, sensor_positions_final)

print("MAP FIM SCORE")
print("FIM det sum:", det_sum)
print("FIM itself:",FIMs)
print(target_localized_successfully, targets, sensor_positions_final, sensor_rad_total, meas_type_total)
new_map = add_targets(ax, targets)


plt.show()
