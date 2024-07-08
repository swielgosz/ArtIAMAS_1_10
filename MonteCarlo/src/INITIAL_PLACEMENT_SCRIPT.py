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
sensor_rad = [50, 50, 50, 50]
sensor_type = ["seismic","acoustic","seismic", "acoustic"]
num_sensors = len(sensor_type)
sensor_comm_ratio = 0.4 # ratio of sensor communication to sensing radius 
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

    sensor_rad, sensor_type, num_sensors, sensor_comm_ratio, meas_type, obj_constraint, step_num = sensor_chars
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
    eig_abs_sum = 0.

    # Set optimizer_var to whatever you want (trace, det, eigenvalue)
    # NOTE: penalties are applied to individual FIMS - see the build_map_FIMS for details!!!!!
    for kk in range(len(FIMs)):
        # FIM correspond to target list one-to-one
        det_mult = det_mult*np.linalg.det(FIMs[kk])
        tr_sum += np.trace(FIMs[kk])
        #eig_abs_sum += np.max(np.linalg.eigvals(FIMs[kk]))

    # Now consider the minimum sensor-sensor distance penalty
    # Need to perform this over all sensors in the WSN
    sens_sens_penalty = min_sensor_distance_penalty(sensor_positions, d_sens_min)
    det_mult = det_mult*sens_sens_penalty #will multiply by zero if sensor distances are not met, 1 if they are
    tr_sum = tr_sum*sens_sens_penalty
    #eig_abs_sum = eig_abs_sum*sens_sens_penalty

    # Now consider the valid placement check
    # Only check over the optimized sensor positions, x. Doesn't make
    # sense to check the pre-existing positions as well
    valid_place_penalty = valid_sensor_penalty(x, terrain)
    det_mult = det_mult*valid_place_penalty #will multiply by zero if sensor distances are not met, 1 if they are
    tr_sum  = tr_sum *valid_place_penalty
    #eig_abs_sum  = eig_abs_sum*valid_place_penalty

    # Finally, check that the sensors are indeed in a WSN and add penalty as necessary
    # Uses _map, constructed above
    valid_WSN = terrain.is_configuration_valid(_map)
    if valid_WSN:
        valid_WSN_penalty = 1 # 1 is satisfaction of WSN
    else:
        # Else, calculate the penalty on WSN
        # Note: only ONE comm radii is accepted. Not written to vary (just yet!)
        # comm_radii should match the inputs into terrain.is_config_valid(_map)
        comm_radii = sensor_rad[0]*sensor_comm_ratio 
        valid_WSN_penalty = WSN_penalty(sensor_positions, comm_radii) 

    det_mult = det_mult*valid_WSN_penalty
    tr_sum = tr_sum*valid_WSN_penalty
    #eig_abs_sum = eig_abs_sum*valid_WSN_penalty

    # For TWO STEP OPTIMIZATION, add an epsilon constraint
    # Now consider the two-step optimization
    def obj_penalty(obj, obj_const):
        mu = 0.25
        # critical point at eps*step1_obj
        return 1/(1+np.exp(-mu*(obj-obj_const)))
    
    if step_num == 1:
        optimizer_var = tr_sum
    elif step_num == 2:
        optimizer_var = det_mult
        # penalize the constraint if not within epsilon
        epsilon = 0.9
        obj_penalty_val = obj_penalty(det_mult, obj_constraint*epsilon)
        #obj_penalty_val = 1
        optimizer_var = optimizer_var*obj_penalty_val

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
target_mesh_pts = [[50, 50],[49, 49],[51, 51], [60, 60], [62, 62], [58, 58], [40, 40], [42, 42], [53, 48], [55, 50], [53, 50], [55, 48]]

for i in range(1):
    step_num = i+1

    if i == 0: # first pass, set constraint to zero
        obj_constraint = 0 

    # Construct tuples to pass in
    sensor_chars = (sensor_rad, sensor_type, num_sensors, sensor_comm_ratio, meas_type, obj_constraint, step_num)
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
        res = optimize.direct(objective_fcn, bounds=bounds, args = sensor_list, vol_tol=vol_tol_test, maxfun = maxfun_1, maxiter = max_iters)
        if i == 0:
            ax_opt.plot(fcn_counter, fcn_eval_list, linestyle='None', marker='+')
        elif i == 1:
            ax_opt.plot(fcn_counter, fcn_eval_list, linestyle='None', marker='o')
        print(res.message, "Final values:", res.x)
        x_out = res.x
        obj_constraint = -1.0*res.fun

    # Multi-objective optimization!
    # If we don't want to run it, just set range(1) instead of (2) in the loop
    if i == 0:
        first_run_positions = x_out
    if i == 1:
        second_run_positions = x_out

    # Save the number of runs we did for later
    num_runs = i+1

# Plot formatting
ax_opt.set_title('Initial Placement Optimization')
ax_opt.set_ylabel('Objective Function Evaluations'), ax_opt.set_xlabel('Iteration Count')
#ax_opt.set_ylim([min(fcn_eval_list)*1.25, max(fcn_eval_list)*1.25])
ax_opt.grid(True, which='minor')  
plt.show()

print('-----------------------------------')

# ---------FORMAT AND ANALYSIS------------------
# Finally, make the map for plotting
# Localize the new map
target_localized_successfully = [1 for _ in range(len(target_mesh_pts ))] #just need a list of 1's 
targets_in = []
for target in target_mesh_pts:
    for pt in target:
        targets_in.append(pt) 

# fcn for computing the trace and det
def return_obj_vals(FIMs):
    det_mult = 1.
    tr_sum = 0.
    for kk in range(len(FIMs)):
        # FIM correspond to target list one-to-one
        det_mult = det_mult*np.linalg.det(FIMs[kk])
        tr_sum += np.trace(FIMs[kk])
    return det_mult, tr_sum

# Analyze first run results
inputs = target_localized_successfully, targets_in, first_run_positions, sensor_rad, meas_type, terrain, 0 #LOS_flag == 1
(FIMs_first_run, det_sum) = build_map_FIMS(inputs)
tr_sum1, det_mult1 = return_obj_vals(FIMs_first_run)
print("First run results (trace and det):", tr_sum1, det_mult1)
print("Associated positions:", first_run_positions)
_map = make_basic_seismic_map(num_sensors, sensor_rad, sensor_type, meas_type, sensor_comm_ratio, first_run_positions)

# Analyze second run results  
if num_runs == 2:
    inputs = target_localized_successfully, targets_in, second_run_positions, sensor_rad, meas_type, terrain, 0 #LOS_flag == 1
    (FIMs_sec_run, det_sum) = build_map_FIMS(inputs)
    tr_sum2, det_mult2 = return_obj_vals(FIMs_sec_run)
    print("First run results (trace and det):", tr_sum2, det_mult2)
    print("Associated positions:", second_run_positions)
    # Reconstructs the map if a second run was made
    _map = make_basic_seismic_map(num_sensors, sensor_rad, sensor_type, meas_type, sensor_comm_ratio, second_run_positions)

# Finally, construct the plot and add the ellipses
ax = terrain.plot_grid(_map)
new_map = add_targets(ax, targets_in)
for i in range(len(targets_in)//2):
    target_i = [targets_in[0+2*i], targets_in[1+2*i]]
    new_map = plot_uncertainty_ellipse(new_map, FIMs_first_run[i], target_i, 2.48, 1, "grey", "analytical")
    if num_runs == 2:
        new_map = plot_uncertainty_ellipse(new_map, FIMs_sec_run[i], target_i, 2.48, 1, "black", "analytical")

# Add some text to the plot to show initial vs final placed sensors
if num_runs == 2:
    for i in range(len(second_run_positions)//2):
        plt.text(second_run_positions[0+2*i]+1, second_run_positions[1+2*i]+1, "In")
elif num_runs == 1:
    for i in range(len(first_run_positions)//2):
        plt.text(first_run_positions[0+2*i]+1, first_run_positions[1+2*i]+1, "In")

plt.show()
