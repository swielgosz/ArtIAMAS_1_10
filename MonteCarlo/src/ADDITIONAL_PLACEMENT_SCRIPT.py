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
from helper_calculations.fisher_information import build_FIM, build_map_FIMS, plot_uncertainty_ellipse, uncertainty_ellipse_stats, return_obj_vals
from helper_calculations.localization_calculations import sensor_localization_routine
from helper_calculations.sensor_vision import get_LOS_coeff
from helper_calculations.penalty_fcn import min_distance_penalty, min_sensor_distance_penalty, valid_sensor_penalty, WSN_penalty
from helper_calculations.target_meshing import make_target_mesh, sample_target_mesh
from helper_calculations.save_vals import save_data

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
### TERRAIN INPUTS ###
terrain_height = 100 #796
terrain_width = 100 #1002
terrain = Configuration(terrain_width, terrain_height)
my_path = Path(__file__).parent / "../data/terrain/GIS_terrain_resize.csv"
terrain.load_from_csv(my_path)

### TARGET INPUTS ###
# Create the target mesh
area_origin =  [[55, 35], [57, 33], [61, 33], [63, 33], [34, 28], [36, 32], [40, 34]] # can add additional origins if using multiple rectangles
area_dim = [[2, 8], [4, 10], [2, 6], [2, 4], [8, 4], [8, 2], [4, 2]] # [width, height]
target_mesh_points = []
step = 1 # more finely sample compared to the initial run!

# Place targets
for i in range(len(area_origin)):
    target_mesh_points = make_target_mesh(target_mesh_points,terrain,area_origin[i], area_dim[i], step)

# Sample from mesh
n = 5 # points to sample from
targets_discovered = sample_target_mesh(target_mesh_points, n)

# Write as a list
targets = []
for target in targets_discovered:
    targets.append(target[0])
    targets.append(target[1])

### SENSOR INPUTS ###
sensor_locs = [54.50633808, 46.32458927, 40.44200985, 39.33369424, 47.27705969, 30.66311726]
sensor_rad = [35, 35, 35]
sensor_type = ["seismic","acoustic","seismic"]
num_sensors = len(sensor_type)
sensor_comm_ratio = 0.75 # ratio of sensor communication to sensing radius 
meas_type = ["radius", "bearing", "radius"]
LOS_flag = 0 # 1 if want to consider LOS, 0 if don't want to

# OPTIMIZATION PARAMETERS
threshold = 4 #minimum distance a sensor must be from a target
    # NOTE: THIS IS ENFORCED IN FISHER_INFORMATION.PY -> BUILD_FIM()!!!
    # GO THERE TO CHANGE THE PARAMETER
d_sens_min = 4 #minimum sensor-sensor distance
maxfun_1 = 50000 #Max fun w/o LOS
max_iters = 50000
printer_counts = 1000 #print the results every ___ fcn calls
# Use list to run a parameter sweep
vol_tol = (1e-5)**(2*2)
# vol_tol = [1e-14] #optimization param
# ------------------------------------------

## DETERMINE THE LOCALIZABLE TARGETS
target_localized_successfully = [1 for _ in targets_discovered]

## NOMINAL FIM MAP SCORE
inputs = target_localized_successfully, targets, sensor_locs, sensor_rad, meas_type, terrain, 0 #LOS_flag == 1
(FIMs_stats, det_sum) = build_map_FIMS(inputs)

for i in range(len(targets)//2):
    target_i = [targets[0+2*i], targets[1+2*i]]
    # Compute the stats
    Cov_Mat, a, b, sx, sy, area = uncertainty_ellipse_stats(FIMs_stats[i], target_i, 1, 1)
    print("Target and area:", target_i, area)

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
    global best_fcn

    counter += 1
    # (0) PRELIMINARIES
    # Initialize the sensor lists
    sensor_rad, sensor_type, meas_type, sensor_positions, num_sensors = [], [], [], [], 0

    # Pull out the old and new sensor lists
    existing_sensor_lists, new_sensor_lists, target_inputs, step_num, obj_constraint = args

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

    # (1) CONSTRUCT THE FIMs
    inputs = target_localized_successfully, targets, sensor_positions, sensor_rad, meas_type, terrain, LOS_flag
    (FIMs, det_sum) = build_map_FIMS(inputs)
    det_sum = 0.
    tr_sum = 0.

    # Set optimizer_var to whatever you want (trace, det, eigenvalue)
    # NOTE: penalties are applied to individual FIMS - see the build_map_FIMS for details!!!!!
    optimizer_var = det_sum
    for kk in range(len(FIMs)):
        # FIM correspond to target list one-to-one
        det_sum = det_sum + np.linalg.det(FIMs[kk])
        tr_sum += np.trace(FIMs[kk])

    # Consider the min target-sensor penalty
    sens_target_penalty = min_distance_penalty(targets, sensor_positions, threshold)
    det_sum = det_sum*(1-sens_target_penalty)
    tr_sum = tr_sum*(1-sens_target_penalty)

    # Now consider the minimum sensor-sensor distance penalty
    # Need to perform this over all sensors in the WSN
    sens_sens_penalty = min_sensor_distance_penalty(sensor_positions, d_sens_min)
    det_sum = det_sum*sens_sens_penalty #will multiply by zero if sensor distances are not met, 1 if they are
    tr_sum = tr_sum*sens_sens_penalty

    # Now consider the valid placement check
    # Only check over the optimized sensor positions, x. Doesn't make
    # sense to check the pre-existing positions as well
    valid_place_penalty = valid_sensor_penalty(x, terrain)
    det_sum = det_sum*valid_place_penalty #will multiply by zero if sensor distances are not met, 1 if they are
    tr_sum = tr_sum*valid_place_penalty

    # Finally, check that the sensors are indeed in a WSN and add penalty as necessary
    # Uses _map, constructed above
    valid_WSN = terrain.is_configuration_valid(_map)
    if valid_WSN:
        valid_WSN_penalty = 1 # 1 is satisfaction of WSN
    else:
        # Else, calculate the penalty on WSN
        # Note: only ONE comm radii is accepted. Not written to vary (just yet!)
        comm_radii = sensor_rad_new[0]*sensor_comm_ratio # just some scalar
        valid_WSN_penalty = WSN_penalty(sensor_positions, comm_radii) 

    det_sum = det_sum*valid_WSN_penalty
    tr_sum = det_sum*valid_WSN_penalty

    # For TWO STEP OPTIMIZATION, add an epsilon constraint
    # Now consider the two-step optimization
    def obj_penalty(obj, obj_const):
        mu = 5
        # critical point at eps*step1_obj
        return 1/(1+np.exp(-mu*(obj-obj_const)))
    
    # PICK OPTIMIZATION ORDER HERE!
    if step_num == 1:
        optimizer_var = tr_sum #CHANGE THIS 
        eps_opt_penalty = 1
        prev_opt = 1
    elif step_num == 2:
        optimizer_var = det_sum # CHANGE THIS
        prev_opt = tr_sum # CHANGE THIS TO MATCH STEP 1
        epsilon = 0.95
        eps_opt_penalty = obj_penalty(prev_opt, epsilon*obj_constraint)

    optimizer_var = optimizer_var*eps_opt_penalty
    
    if optimizer_var > best_fcn:
        best_fcn = optimizer_var.copy()

    # Maximize!!
    if counter % printer_counts == 0:
        print(counter, optimizer_var, det_sum, tr_sum, eps_opt_penalty)

    fcn_eval_list.append(best_fcn)
    fcn_counter.append(counter)
    return -optimizer_var # Minimize the negative det(FIMS)

# ------------------------------------------
# Set a list of additional sensors to place down
# SENSOR LIST
sensor_rad_new = [35, 35]
sensor_type_new = ["seismic", "acoustic"]
num_sensors_new = len(sensor_type_new)
sensor_comm_ratio_new = sensor_comm_ratio # ratio of sensor communication to sensing radius 
meas_type_new = ["radius", "radius"]

# RUN THE OPTIMIZER!
optimized_vals = []

for i in range(2):
    maxfun = maxfun_1
    step_num = i+1
    
    obj_constraint = 0

    # Construct tuples to pass in
    LOS_flag = 0
    existing_sensor_lists = (sensor_locs, sensor_rad, sensor_type, num_sensors, sensor_comm_ratio, meas_type)
    new_sensor_lists = (sensor_rad_new, sensor_type_new, num_sensors_new, sensor_comm_ratio_new, meas_type_new)
    target_ins = (targets, target_localized_successfully, terrain)

    # Set the bounds (controls the number of variables to reason over)
    bounds = [(0, terrain_width-1),(0, terrain_height-1)]*(num_sensors_new) #bounds placing points to only on the map
    sensor_list = (existing_sensor_lists, new_sensor_lists, target_ins, step_num, obj_constraint)   

    # Call the optimizer!
    counter, best_fcn = 0, 0

    # This part of the script just plots the obj fcn vs count #
    # so we can check progress
    fig = plt.figure()
    ax_opt = fig.add_subplot()
    fcn_eval_list, fcn_counter = [], []
    print("Start optimizer")
    res = optimize.direct(objective_fcn, bounds=bounds, args = sensor_list, vol_tol=vol_tol, maxfun = maxfun, maxiter = max_iters)
    ax_opt.plot(fcn_counter, fcn_eval_list, linestyle='None', marker='+')
    print(res.message, "Final values:", res.x)
    x_out = res.x

    res = optimize.minimize(objective_fcn, x_out, args = sensor_list, method='BFGS', jac='3-point', options={'gtol': 0.01})
    obj_constraint = -1.0*res.fun
    print("final objective val:", res.fun)

    optimized_vals.append(x_out)
    ax_opt.set_title('Function Evaluation across Optimization - No LOS considerations')


print('-----------------------------------')
print("Final positions:", x_out)

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
for j, sensor_positions in enumerate(optimized_vals):

    sensors_in = sensor_locs.copy()
    for val in sensor_positions:
        sensors_in.append(val)

    print("Sensor list:", sensors_in)
    inputs = target_localized_successfully, targets, sensors_in, sensor_rad_total, meas_type_total, terrain, 0 #LOS_flag == 0
    (FIMs_final, det_sum) = build_map_FIMS(inputs)
    det_mult1, tr_sum1 = return_obj_vals(FIMs_final)
    print(i+1," Placement det and trace sums:", det_sum, tr_sum1)
    for i in range(len(targets)//2):
        target_i = [targets[0+2*i], targets[1+2*i]]
        # Compute the stats
        Cov_Mat, a, b, sx, sy, area = uncertainty_ellipse_stats(FIMs_final[i], target_i, 1, 1)
        print("Target and area:", target_i, area)


# Finally, construct the plot and add the ellipses
_map = make_basic_seismic_map(num_sens_total, sensor_rad_total, sensor_type_total, meas_type_total, sensor_comm_ratio, sensors_in)
ax = terrain.plot_grid(_map)
new_map = add_targets(ax, targets)
color = ["blue", "black"]
for j, sensor_positions in enumerate(optimized_vals):
    # Generate the sensor list
    sensors_in = sensor_locs.copy()
    for val in sensor_positions:
        sensors_in.append(val)

    for i in range(len(targets)//2):
        target_i = [targets[0+2*i], targets[1+2*i]]
        new_map = plot_uncertainty_ellipse(new_map, FIMs_final[i], target_i, 1, 1, color[j], "analytical")
        #new_map = plot_uncertainty_ellipse(new_map, FIMs_LOS[i], target_i, 2.48, 1, "black", "analytical")

# Plot the initial uncert ellipses
for i in range(len(targets)//2):
        target_i = [targets[0+2*i], targets[1+2*i]]
        new_map = plot_uncertainty_ellipse(new_map, FIMs_stats[i], target_i, 1, 1, "grey", "analytical")

# Add some text to the plot to show initial vs final placed sensors
for i in range(len(sensor_locs)//2):
    plt.text(sensor_locs[0+2*i]+1, sensor_locs[1+2*i]+1, "In")

for i in range(len(x_out)//2):
    plt.text(x_out[0+2*i]+1, x_out[1+2*i]+1, "F") #plots the LOS sensors

# save off sensor information
sensor_save_off = [sensors_in, meas_type_total]
sensor_save_names = ["sensor locations", "sensor_types"]
save_data(sensor_save_off, sensor_save_names, "final_sensor_placements_Case1", "Case 1 final sensor placement")

# save off FIMs
FIM_save_off = [targets, FIMs_stats, FIMs_final] #targets, initial fims, final fims
FIM_save_name = ["targets sampled", "initial FIMs", "final FIMs"]
save_data(FIM_save_off, FIM_save_name, "final_sensor_FIMs_Case1", "Case 1 final FIMs")


plt.show()
