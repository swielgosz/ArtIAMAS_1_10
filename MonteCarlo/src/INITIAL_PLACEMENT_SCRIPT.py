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
# from helper_calculations.sensor_vision import sensor_vision, sensor_reading, target_localization_both, target_localization_bearing, target_localization_radius
from helper_calculations.fisher_information import build_FIM, build_map_FIMS, plot_uncertainty_ellipse, return_obj_vals, uncertainty_ellipse_stats
from helper_calculations.localization_calculations import sensor_localization_routine
from helper_calculations.sensor_vision import get_LOS_coeff
from helper_calculations.penalty_fcn import min_distance_penalty, min_sensor_distance_penalty, valid_sensor_penalty, WSN_penalty
from helper_calculations.save_vals import save_data
from helper_calculations.target_meshing import make_target_mesh

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
sensor_comm_ratio = 0.5 # ratio of sensor communication to sensing radius 
meas_type = ["radius", "bearing", "radius", "bearing"] #radius or bearing
LOS_flag = 0 # 1 if want to consider LOS, 0 if don't want to

# OPTIMIZATION PARAMETERS
threshold = 6 #minimum distance a sensor must be from a target
d_sens_min = 4 #minimum sensor-sensor distance
printer_counts = 50 #print the results every ___ fcn calls

# DIRECT Parameters
vol_tol = 1e-30
max_iters = 15000
maxfun_1 = 15000

# DE Parameters - optimized!
popsize = 7
mutation = 0.8
recombination = 0.8

# Sim Annealing Params
initial_temp = 3000
restart_temp_ratio = 0.0005
visit = 2.75

# ------------------------------------------

# ----------INITIAL PLACEMENT---------------
# Now define the objective fcn
def objective_fcn(x, *args):
    global counter
    global fcn_eval_list
    global fcn_counter
    global best_fcn

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

    # Compute the FIMs
    inputs = target_localized_successfully, targets_in, sensor_positions, sensor_rad, meas_type, terrain, LOS_flag
    (FIMs, det_sum) = build_map_FIMS(inputs)
    det_mult = 1.
    tr_sum = 0.
    eig_abs_sum = 0.

    # Set optimizer_var to whatever you want (trace, det, eigenvalue)
    # NOTE: penalties are applied to individual FIMS - see the build_map_FIMS for details!!!!!
    for kk in range(len(FIMs)):
        # FIM correspond to target list one-to-one
        det_mult = det_mult + np.linalg.det(FIMs[kk])
        tr_sum += np.trace(FIMs[kk])
        #eig_abs_sum += np.max(np.linalg.eigvals(FIMs[kk]))

    # Consider the minimum sensor-target penalty
    sens_target_penalty = min_distance_penalty(targets_in, sensor_positions, threshold)
    det_mult = det_mult*sens_target_penalty #will multiply by zero if sensor distances are not met, 1 if they are
    tr_sum = tr_sum*sens_target_penalty

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
    tr_sum  = tr_sum*valid_place_penalty
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
        mu = 1
        # critical point at eps*step1_obj
        return 1/(1+np.exp(-mu*(obj-obj_const)))
    
    # PICK OPTIMIZATION ORDER HERE!
    if step_num == 1:
        optimizer_var = tr_sum #CHANGE THIS 
        eps_opt_penalty = 1
        prev_opt = 1
    elif step_num == 2:
        optimizer_var = det_mult # CHANGE THIS
        prev_opt = tr_sum # CHANGE THIS TO MATCH STEP 1
        epsilon = 0.95
        eps_opt_penalty = obj_penalty(prev_opt, epsilon*obj_constraint)

    optimizer_var = optimizer_var*eps_opt_penalty
    
    #optimizer_var = optimizer_var * eps_opt_penalty
    if optimizer_var > best_fcn:
        best_fcn = optimizer_var.copy()

    # If a valid config, return (-1)det(FIMs)
    if valid_placement_check:# and valid_WSN:
        # Maximize the determinant of the map
        if counter % printer_counts == 0:
            print(counter, tr_sum, det_mult, sens_sens_penalty, valid_WSN_penalty, sens_target_penalty, valid_place_penalty)

        fcn_eval_list.append(best_fcn)
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
# Initialize to save off data later as a csv
data_save_off = []
data_name_list = []

# Generate target mesh here!
#target_mesh_pts = [[50, 50],[49, 49],[51, 51], [60, 60], [62, 62], [58, 58], [40, 40], [42, 42], [53, 48], [55, 50], [53, 50], [55, 48]]

# Define rectangular areas to place targets in
num_areas = 1
area_origin = [[44, 42]] # can add additional origins if using multiple rectangles
area_dim = [[6,6]] # same as above
target_mesh_points = []

# Place targets
for i in range(num_areas):
    target_mesh_pts = make_target_mesh(target_mesh_points,terrain,area_origin[i-1], area_dim[i-1])

# Define the list of optimizers being tested
# options are (exactly): 'DIRECT', 'DE', 'SA'
optimizers = ['DIRECT']#, 'DE', 'SA']

optimized_vals = []

# Run the multi-objective placements!
# To pick what's being optimized first, change it around in the objective fcn
for i in range(1): # set to one if doing single-step. Two otherwise
    # Generate figs and plots
    fig = plt.figure()
    ax_opt = fig.add_subplot()
    step_num = i+1

    if step_num == 2: #overwrite the prev optimized vals if two-step
        optimized_vals = []

    if i == 0: # first pass, set constraint to zero
        obj_constraint = 0 

    # Construct tuples to pass in
    sensor_chars = (sensor_rad, sensor_type, num_sensors, sensor_comm_ratio, meas_type, obj_constraint, step_num)
    target_ins = (target_mesh_pts, terrain)

    # Set the bounds (controls the number of variables to reason over)
    bounds = [(0, terrain_width-1),(0, terrain_height-1)]*(num_sensors) #bounds placing points to only on the map
    sensor_list = (sensor_chars, target_ins)   

    # Loop over the tested optimizers
    for optimizer in optimizers: #lets you run multiple tolerances if you want
        fcn_eval_list, fcn_counter = [], []
        counter, best_fcn = 0, 0
        
        print("Start", optimizer, "optimizer")

        if optimizer == 'DIRECT':
            res = optimize.direct(objective_fcn, bounds=bounds, args = sensor_list, vol_tol=vol_tol, maxfun = maxfun_1, maxiter = max_iters)
            data_name_list.append("DIRECT Optimizer curve")
        elif optimizer == 'DE':
            res = optimize.differential_evolution(objective_fcn, bounds=bounds, args=sensor_list, strategy='best1bin', 
                                                  popsize=popsize, tol=0.01, mutation=mutation, recombination=recombination, maxiter=(int(maxfun_1/popsize/num_sensors)-1))
            data_name_list.append("Diff Ev Optimizer curve")
        elif optimizer == 'SA':
            res = optimize.dual_annealing(objective_fcn, bounds=bounds, args=sensor_list, maxfun=maxfun_1, maxiter=max_iters, no_local_search=True, 
                                          initial_temp=initial_temp, restart_temp_ratio=restart_temp_ratio, visit=visit, accept = -0.01)
            data_name_list.append("Sim Anneal Optimizer curve")
        
        print("Complete optimizing under", optimizer, "with value and iterations", res.fun, counter)
        print(res.message)
        #print("First optimizer final values:", res.x)
        
        # Now try fine tuning by optimizing locally
        x_out = res.x #save off prev results to start fine-tuning
        res = optimize.minimize(objective_fcn, x_out, args = sensor_list, method='BFGS', jac='3-point', options={'gtol': 0.005})
        obj_constraint = -1.0*res.fun
        print("Complete fine optimizing under", optimizer, "with value and iterations", res.fun, counter)

        # Format the plots
        ax_opt.plot(fcn_counter, fcn_eval_list, linestyle='None', marker='+', label = optimizer)

        # Save off optimized vals
        optimized_vals.append(res.x)
        x_out = res.x
        
        data_save_off.append(fcn_eval_list)
        data_save_off.append(fcn_counter)
        data_name_list.append("steps in solve")
    
    # End for optimizer ... loop

    # Save the number of runs we did for later
    num_runs = i+1

# Save data off

# Plot formatting
ax_opt.set_ylabel('Objective Function Evaluations'), ax_opt.set_xlabel('Iteration Count')
#ax_opt.set_ylim([min(fcn_eval_list)*1.25, max(fcn_eval_list)*1.25])
ax_opt.grid(True, which='minor')  
plt.legend()
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

# Analyze results
for i, sensor_list in enumerate(optimized_vals):
    inputs = target_localized_successfully, targets_in, sensor_list, sensor_rad, meas_type, terrain, 0 #LOS_flag == 1
    (FIMs_first_run, det_sum) = build_map_FIMS(inputs)
    det_mult1, tr_sum1 = return_obj_vals(FIMs_first_run)
    print(optimizers[i], " run results (trace and det):", tr_sum1, det_mult1)

print('-----------------------------------')
# ---------PLOTTING AND VISUALIZATION------------------
# Comment out the random ones if we don't want to visualize them! 
# It's not always necessary because the ellipses can get unecessarily large

sensor_save_data, sensor_save_names = [], []
FIM_save_data, FIM_save_names = [], []

_map = make_basic_seismic_map(num_sensors, sensor_rad, sensor_type, meas_type, sensor_comm_ratio, x_out)
ax = terrain.plot_grid(_map)
new_map = add_targets(ax, targets_in)

# Set colors for optimizer. Must be the same length as number of optimizers
colors = ["grey", "blue", "red"]

FIM_save_data.append(targets_in)
FIM_save_names.append("targets calculated")

for j, optimizer in enumerate(optimizers):
    sensor_list = optimized_vals[j]
    inputs = target_localized_successfully, targets_in, sensor_list, sensor_rad, meas_type, terrain, 0 #LOS_flag == 1
    (FIMs_stats, det_sum) = build_map_FIMS(inputs)
    FIM_areas, maj_axis= [], []

    for i in range(len(targets_in)//2):
        target_i = [targets_in[0+2*i], targets_in[1+2*i]]
        new_map = plot_uncertainty_ellipse(new_map, FIMs_stats[i], target_i, 1, 1, colors[j], "analytical")
        #new_map = plot_uncertainty_ellipse(new_map, FIMs_random[i], target_i, 2.48, 1, "blue", "numerical")
        
        # Compute the stats
        Cov_Mat, a, b, sx, sy, area_first = uncertainty_ellipse_stats(FIMs_stats[i], target, 2.48, 1)
        FIM_areas.append(area_first)
        maj_axis.append(a)

    print(optimizer, " mean area:", np.mean(FIM_areas))
    print(optimizer, " largest axis", np.max(maj_axis))

    # Save off data!
    sensor_save_data.append(sensor_list)
    sensor_save_names.append(optimizer)

    FIM_save_data.append(FIM_areas)
    FIM_save_data.append(maj_axis)
    FIM_save_names.append(optimizer+" mean area")
    FIM_save_names.append(optimizer+" largest axis")


# Add to plot
for j, optimizer in enumerate(optimizers):
    print(optimizer, optimized_vals[j])
    for i in range(len(optimized_vals[j])//2):
            ax.plot(optimized_vals[j][0+2*i], optimized_vals[j][1+2*i], marker='x', color = colors[j])
            ax.text(optimized_vals[j][0+2*i]+1, optimized_vals[j][1+2*i]+1, optimizer)

# Save off data!
save_data(sensor_save_data, sensor_save_names, "init_sensor_placements", "Initial sensor placement results for differnet optimizers")
save_data(FIM_save_data, FIM_save_names, "init_fim_characteristics", "Initial sensor placement FIM characteristics")
save_data(data_save_off, data_name_list, "Init_place_opt_evals", "Initial Placement Optimization Runs")

plt.show()
