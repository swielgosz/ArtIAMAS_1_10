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
from helper_calculations.penalty_fcn import min_distance_penalty, min_sensor_distance_penalty

# -------- DESCRIPTION -----------------
# This script was makes an ADDITIONAL of a set of sensors by minimizing
# the determinant of a FIM constructed about the entire map. Optimization is
# performed numerically using the DIviding RECTangles algorithm with constraints
# Define the initial sensors in the PARAMETERS section
# Then define the additional sensors after the objective function section
# The script will output a plot with the sensors and the before and after 
# uncertainty ellipses
# --------------------------------------

# Add in a penalty function!!
# JM edits to show branch functionality

# --------------PARAMETERS------------------
# TERRAIN INPUTS
terrain_height = 100 #796
terrain_width = 100 #1002
terrain = Configuration(terrain_width, terrain_height)
my_path = Path(__file__).parent / "../data/terrain/GIS_terrain_resize.csv"
terrain.load_from_csv(my_path)

# INITIAL SENSOR LIST
# Characteristic list
sensor_rad = [25, 25, 25]
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

sens1 = [29, 44] 
sens2 = [33, 29] 
sens3 = [58, 34] 
#sens4 = [44, 29] 
sensor_locs.extend(sens1)
sensor_locs.extend(sens2)
sensor_locs.extend(sens3)
#sensor_locs.extend(sens4)


# OPTIMIZATION PARAMETERS
threshold = 6 #minimum distance a sensor must be from a target
d_sens_min = 3 #minimum sensor-sensor distance
vol_tol = 1e-18
maxfun_1 = 50000 #Max fun w/o LOS
maxfun_2 = 50000 #Max fun w/ LOS
max_iters = 5000
printer_counts = 1000 #print the results every ___ fcn calls
# Use list to run a parameter sweep
vol_tol = [1e-30]
# vol_tol = [1e-14] #optimization param
# ------------------------------------------

## DETERMINE THE LOCALIZABLE TARGETS
inputs = sensor_rad, meas_type, sensor_locs, targets, None
(target_localized_successfully, _) = sensor_localization_routine(inputs)
target_localized_successfully = [1,1,1,1]

## NOMINAL FIM MAP SCORE
inputs = target_localized_successfully, targets, sensor_locs, sensor_rad,meas_type, terrain, LOS_flag
(FIMs, det_sum) = build_map_FIMS(inputs)
print(FIMs)
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
    targets, target_localized_successfully, terrain = target_inputs

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
        i, j = x[0+2*k], x[1+2*k] #let these vary and not be ints
        sensor_positions.extend([i, j])

    # Create a map with the existing sensors and new sensor
    _map = make_basic_seismic_map(num_sensors, sensor_rad, sensor_type, meas_type, sensor_comm_ratio, sensor_positions)

    # (1) CHECK FOR VALID CONFIGURATIONS
    valid_placement_check, valid_WSN = True, True

    # (1.1) First check for valid placement on additionally placed sensors
    #for check in range(len(x)//2):
    #    i, j = int(x[0+2*check]), int(x[1+2*check]) #convert to ints 

        # If ANY sensors are invalid, the entire config is invalid
        # Note that (i,j) are switched to respect dims of terrain class
        # this script is (x,y) = width, height. 
        # valid sensor check is (x, y) = height, width
    #    if not terrain.is_valid_sensor_location(j, i):
    #        valid_placement_check = False
    
    # (1.2) Now check for WSN
    #valid_WSN = terrain.is_configuration_valid(_map)

    # (1.3) Ensure the bearing targets are a minimum distance away
    # Otherwise the 1/r^2 will cause it to just blow up
    # Sweep over all the additionally-placed sensors
    #for check in range(len(x)//2):
    #    i, j = int(x[0+2*check]), int(x[1+2*check])

        # Compute the distance from it to a target
    #    for k in range(len(targets)//2):
    #        tx, ty = targets[0+2*k], targets[1+2*k]
    #        dist = ((tx-i)**2+(ty-j)**2)**(1/2)
            # If w/in the limit, this placement is invalid!
    #        if dist < threshold:
    #            valid_placement_check = False

    # (2) CALCULATE THE SCORE
    # (2.2) Construct the FIM's per target + calculate the det score
    #sensor_positions[-2] = x[0]
    #sensor_positions[-1] = x[1]

    inputs = target_localized_successfully, targets, sensor_positions, sensor_rad, meas_type, terrain, LOS_flag
    (FIMs, det_sum) = build_map_FIMS(inputs)
    det_mult = 1.
    tr_sum = 0.

    # Set optimizer_var to whatever you want (trace, det, eigenvalue)
    # NOTE: penalties are applied to individual FIMS - see the build_map_FIMS for details!!!!!
    optimizer_var = det_mult
    for kk in range(len(FIMs)):
        # FIM correspond to target list one-to-one
        optimizer_var = optimizer_var*np.linalg.det(FIMs[kk])
        #optimizer_var += np.trace(FIMs[kk])

    # Now consider the minimum sensor-sensor distance penalty
    sens_sens_penalty = min_sensor_distance_penalty(sensor_positions, d_sens_min)
    optimizer_var = optimizer_var*sens_sens_penalty #will multiply by zero if sensor distances are not met, 1 if they are

    # If a valid config, return (-1)det(FIMs)
    if valid_placement_check and valid_WSN:
        # Maximize the determinant of the map
        if counter % printer_counts == 0:
            print(counter, optimizer_var, sensor_positions, FIMs)

        fcn_eval_list.append(optimizer_var)
        fcn_counter.append(counter)
        return -optimizer_var # Minimize the negative det(FIMS)
    
    # If an invalid construction, then fail and return 0
    else:
        fcn_eval_list.append(0)
        fcn_counter.append(counter)
        if counter % printer_counts == 0:
            print(counter, 0, sensor_positions, valid_WSN, valid_placement_check)
        return 0

# ------------------------------------------
# Set a list of additional sensors to place down
# SENSOR LIST
sensor_rad_new = [25, 25, 25]
sensor_type_new = ["seismic", "seismic", "acoustic"]
num_sensors_new = len(sensor_type_new)
sensor_comm_ratio_new = sensor_comm_ratio # ratio of sensor communication to sensing radius 
meas_type_new = ["radial", "radial","bearing"]

# Now we run the optimizer w/ and w/o LOS considerations
final_pos_LOS = []
final_pos_no_LOS = []

for i in range(2):
    if i == 0:
        maxfun = maxfun_1
    elif i == 1:
        maxfun = maxfun_2
    # Construct tuples to pass in
    # CHANGE BACK WHEN DONE!!!!
    LOS_flag = 0
    existing_sensor_lists = (sensor_locs, sensor_rad, sensor_type, num_sensors, sensor_comm_ratio, meas_type)
    new_sensor_lists = (sensor_rad_new, sensor_type_new, num_sensors_new, sensor_comm_ratio_new, meas_type_new)
    target_ins = (targets, target_localized_successfully, terrain)

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

    if i == 0:
        final_pos_no_LOS = sensor_locs.copy()
        final_pos_no_LOS.extend(x_out)
        ax_opt.set_title('Function Evaluation across Optimization - No LOS considerations')
    if i == 1:
        final_pos_LOS = sensor_locs.copy()
        final_pos_LOS.extend(x_out)
        ax_opt.set_title('Function Evaluation across Optimization - LOS considerations')

print('-----------------------------------')
print("w/o LOS positions:", final_pos_no_LOS)
print("w/ LOS positions :", final_pos_LOS)

# ------------------------------------------
# Formatting
ax_opt.set_ylabel('Objective Function Evaluations'), ax_opt.set_xlabel('Iteration Count')
#ax_opt.set_ylim([min(fcn_eval_list)*1.25, max(fcn_eval_list)*1.25])
ax_opt.grid(True, which='minor')  
plt.show()

# Create the finalized lists and create a new map
num_sens_total = num_sensors+num_sensors_new
sensor_rad_total, sensor_type_total, meas_type_total = [], [], []
for i in range(num_sensors):
    sensor_rad_total.append(sensor_rad[i])
    sensor_type_total.append(sensor_type[i])
    meas_type_total.append(meas_type[i])

for i in range(num_sensors_new):
    sensor_rad_total.append(sensor_rad_new[i])
    sensor_type_total.append(sensor_type_new[i])
    meas_type_total.append(meas_type_new[i])


# Finally, make the map for plotting
# Localize the new map
inputs = target_localized_successfully, targets, final_pos_no_LOS, sensor_rad_total, meas_type_total, terrain, 1 #LOS_flag == 1
(FIMs_no_LOS, det_sum) = build_map_FIMS(inputs)
print("No LOS Considerations in planning:", det_sum)
print("No LOS Considerations in planning:", FIMs_no_LOS)

inputs = target_localized_successfully, targets, final_pos_LOS, sensor_rad_total, meas_type_total, terrain, 1
(FIMs_LOS, det_sum) = build_map_FIMS(inputs)
print("LOS Considerations in planning:", det_sum)
print("LOS Considerations in planning:", FIMs_LOS)

# Finally, construct the plot and add the ellipses
_map = make_basic_seismic_map(num_sens_total, sensor_rad_total, sensor_type_total, meas_type_total, sensor_comm_ratio, final_pos_LOS)
ax = terrain.plot_grid(_map)
new_map = add_targets(ax, targets)
for i in range(len(targets)//2):
    target_i = [targets[0+2*i], targets[1+2*i]]
    new_map = plot_uncertainty_ellipse(new_map, FIMs_no_LOS[i], target_i, 2.48, 1, "grey", "analytical")
    new_map = plot_uncertainty_ellipse(new_map, FIMs_LOS[i], target_i, 2.48, 1, "black", "analytical")

# Add some text to the plot to show initial vs final placed sensors
for i in range(len(sensor_locs)//2):
    plt.text(sensor_locs[0+2*i]+1, sensor_locs[1+2*i]+1, "In")

for i in range(len(x_out)//2):
    plt.text(x_out[0+2*i]+1, x_out[1+2*i]+1, "F") #plots the LOS sensors

plt.show()
