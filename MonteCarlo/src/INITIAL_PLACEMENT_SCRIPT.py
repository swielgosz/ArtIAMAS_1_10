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
from helper_calculations.penalty_fcn import min_distance_penalty, min_sensor_distance_penalty, valid_sensor_penalty, WSN_penalty

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
sensor_rad = [100, 100, 100, 100]
sensor_type = ["seismic","acoustic","seismic", "acoustic"]
num_sensors = len(sensor_type)
sensor_comm_ratio = 1.5 # ratio of sensor communication to sensing radius 
meas_type = ["radius", "bearing", "radius", "bearing"] #radius or bearing
LOS_flag = 0 # 1 if want to consider LOS, 0 if don't want to

# OPTIMIZATION PARAMETERS
threshold = 8 #minimum distance a sensor must be from a target
    # NOTE: THIS IS ENFORCED IN FISHER_INFORMATION.PY -> BUILD_FIM()!!!
    # GO THERE TO CHANGE THE PARAMETER
d_sens_min = 4 #minimum sensor-sensor distance
vol_tol = 1e-18
maxfun_1 = 50000 #Max fun w/o LOS
max_iters = 20000
printer_counts = 5000 #print the results every ___ fcn calls
# Use list to run a parameter sweep
vol_tol = [1e-30]
# vol_tol = [1e-14] #optimization param
# ------------------------------------------

# ----------INITIAL PLACEMENT---------------
# Now define the objective fcn
def objective_fcn(x, *args):
    global counter
    global fcn_eval_list
    global fcn_counter

    counter += 1
    # (0) PRELIMINARIES
    # Pull out the old and new sensor lists
    sensor_chars, target_inputs = args

    sensor_rad, sensor_type, num_sensors, sensor_comm_ratio, meas_type = sensor_chars
    target_mesh, terrain = target_inputs
    
    sensor_positions = []
    # Add the optimization variables to the lists
    for k in range(len(x)//2):
        i, j = x[0+2*k], x[1+2*k] #let these vary and not be ints
        sensor_positions.extend([i, j])

    # Create a map with the existing sensors and new sensor
    _map = make_basic_seismic_map(num_sensors, sensor_rad, sensor_type, meas_type, sensor_comm_ratio, sensor_positions)

    # (1) CHECK FOR VALID CONFIGURATIONS
    # Relic of the past that I'm too lazy to change - just leave as true
    valid_placement_check, valid_WSN = True, True

    target_localized_successfully = [1 for _ in range(len(target_mesh))] #just need a list of 1's
    
    # Break out the mesh pts:
    targets_in = []
    for target in target_mesh:
        for pt in target:
            targets_in.append(pt) 

    inputs = target_localized_successfully, targets_in, sensor_positions, sensor_rad, meas_type, terrain, LOS_flag
    (FIMs, det_sum) = build_map_FIMS(inputs)
    det_mult = 1.
    tr_sum = 0.

    # Set optimizer_var to whatever you want (trace, det, eigenvalue)
    # NOTE: penalties are applied to individual FIMS - see the build_map_FIMS for details!!!!!
    optimizer_var = tr_sum
    for kk in range(len(FIMs)):
        # FIM correspond to target list one-to-one
        optimizer_var = optimizer_var*np.linalg.det(FIMs[kk])
        #optimizer_var += np.trace(FIMs[kk])

    # Now consider the minimum sensor-sensor distance penalty
    # Need to perform this over all sensors in the WSN
    sens_sens_penalty = min_sensor_distance_penalty(sensor_positions, d_sens_min)
    optimizer_var = optimizer_var*sens_sens_penalty #will multiply by zero if sensor distances are not met, 1 if they are

    # Now consider the valid placement check
    # Only check over the optimized sensor positions, x. Doesn't make
    # sense to check the pre-existing positions as well
    valid_place_penalty = valid_sensor_penalty(x, terrain)
    optimizer_var = optimizer_var*valid_place_penalty #will multiply by zero if sensor distances are not met, 1 if they are

    # Finally, check that the sensors are indeed in a WSN and add penalty as necessary
    # Uses _map, constructed above
    valid_WSN = terrain.is_configuration_valid(_map)
    if valid_WSN:
        valid_WSN_penalty = 1 # 1 is satisfaction of WSN
    else:
        # Else, calculate the penalty on WSN
        # Note: only ONE comm radii is accepted. Not written to vary (just yet!)
        comm_radii = sensor_rad[0]*sensor_comm_ratio # just some scalar
        valid_WSN_penalty = WSN_penalty(sensor_positions, comm_radii) 

    optimizer_var = optimizer_var*valid_WSN_penalty

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

# ---------RUN THE OPTIMIZER-------------------

# Generate target mesh here!
target_mesh_pts = [[50, 50],[49, 49],[51, 51], [60, 60], [62, 62], [58, 58], [40, 40], [42, 42]]

for i in range(1):
    maxfun = maxfun_1

    # Construct tuples to pass in
    LOS_flag = 0
    
    sensor_chars = (sensor_rad, sensor_type, num_sensors, sensor_comm_ratio, meas_type)
    target_ins = (target_mesh_pts, terrain)

    # Set the bounds (controls the number of variables to reason over)
    bounds = [(0, terrain_width-1),(0, terrain_height-1)]*(num_sensors) #bounds placing points to only on the map
    sensor_list = (sensor_chars, target_ins)   

    # Call the optimizer!
    counter = 0

    # This part of the script just plots the obj fcn vs count #
    # so we can check progress
    fig = plt.figure()
    ax_opt = fig.add_subplot()
    for count, vol_tol_test in enumerate(vol_tol): #lets you run multiple tolerances if you want
        fcn_eval_list, fcn_counter = [], []
        print("Start optimizer")
        res = optimize.direct(objective_fcn, bounds=bounds, args = sensor_list, vol_tol=vol_tol_test, maxfun = maxfun, maxiter = max_iters)
        ax_opt.plot(fcn_counter, fcn_eval_list, linestyle='None', marker='+')
        print(res.message, "Final values:", res.x)
        x_out = res.x

    # Plot formatting
    ax_opt.set_title('Initial Placement Optimization')
    ax_opt.set_ylabel('Objective Function Evaluations'), ax_opt.set_xlabel('Iteration Count')
    #ax_opt.set_ylim([min(fcn_eval_list)*1.25, max(fcn_eval_list)*1.25])
    ax_opt.grid(True, which='minor')  
    plt.show()
    
    # Final optimized values!
    final_sensor_locs = x_out

print('-----------------------------------')

# ---------FORMAT AND ANALYSIS------------------
# Formatting


# Finally, make the map for plotting
# Localize the new map
target_localized_successfully = [1 for _ in range(len(target_mesh_pts ))] #just need a list of 1's 
targets_in = []
for target in target_mesh_pts:
    for pt in target:
        targets_in.append(pt) 

inputs = target_localized_successfully, targets_in, final_sensor_locs, sensor_rad, meas_type, terrain, 0 #LOS_flag == 1
(FIMs_no_LOS, det_sum) = build_map_FIMS(inputs)
print("No LOS Considerations in planning:", det_sum)
print("No LOS Considerations in planning:", FIMs_no_LOS)

#inputs = target_localized_successfully, targets, final_pos_LOS, sensor_rad_total, meas_type_total, terrain, 1
#(FIMs_LOS, det_sum) = build_map_FIMS(inputs)
#print("LOS Considerations in planning:", det_sum)
#print("LOS Considerations in planning:", FIMs_LOS)

# Finally, construct the plot and add the ellipses
_map = make_basic_seismic_map(num_sensors, sensor_rad, sensor_type, meas_type, sensor_comm_ratio, final_sensor_locs)
ax = terrain.plot_grid(_map)
new_map = add_targets(ax, targets_in)
for i in range(len(targets_in)//2):
    target_i = [targets_in[0+2*i], targets_in[1+2*i]]
    new_map = plot_uncertainty_ellipse(new_map, FIMs_no_LOS[i], target_i, 2.48, 1, "grey", "analytical")
    #new_map = plot_uncertainty_ellipse(new_map, FIMs_LOS[i], target_i, 2.48, 1, "black", "analytical")

# Add some text to the plot to show initial vs final placed sensors
for i in range(len(final_sensor_locs)//2):
    plt.text(final_sensor_locs[0+2*i]+1, final_sensor_locs[1+2*i]+1, "In")

plt.show()
