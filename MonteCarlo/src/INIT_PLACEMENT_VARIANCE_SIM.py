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
from helper_calculations.target_meshing import make_target_mesh, upper_achieveable_limit

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
sensor_type = ["seismic","seismic","acoustic","acoustic"]
num_sensors = len(sensor_type)
sensor_comm_ratio = 0.5 # ratio of sensor communication to sensing radius 
meas_type = ["radius", "radius", "bearing", 'bearing']
LOS_flag = 0 # 1 if want to consider LOS, 0 if don't want to

# OPTIMIZATION PARAMETERS
threshold = 4 #minimum distance a sensor must be from a target
d_sens_min = 4 #minimum sensor-sensor distance
printer_counts = 5000 #print the results every ___ fcn calls

# DIRECT Parameters
vol_tol = 1e-36 # should scale (about) as (1e-4)^(2N) where N is the num of sensors
max_iters = 25000
maxfun_1 = 20000

# DE Parameters - optimized!
popsize = 7
mutation = 0.8
recombination = 0.8
max_fun_DE = 15000

# Sim Annealing Params
initial_temp = 3000
restart_temp_ratio = 0.0005
visit = 2.75
max_iter_SA = 12000

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
    det_mult = det_mult*(1-sens_target_penalty) #will multiply by zero if sensor distances are not met, 1 if they are
    tr_sum = tr_sum*(1-sens_target_penalty)

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
    #valid_WSN = terrain.is_configuration_valid(_map)
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
        optimizer_var = det_mult #CHANGE THIS 
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
            print(counter, best_fcn, tr_sum, det_mult, sens_sens_penalty, valid_WSN_penalty, sens_target_penalty, valid_place_penalty, eps_opt_penalty)

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

# Areas
area_list = [10] 
area_plot = [l*l for l in area_list]

# Data to compare
det_area_list, tr_area_list = [], []
det_largest_axis, tr_largest_axis = [], []
det_mean_var, tr_mean_var = [], []
det_vals, tr_vals = [], []
opt_var = []
sensor_positions_save = []
sensor_names_save = []

for area_size in area_list: # For EACH tested area

    # Construct the mesh
    area_dim, area_origin = [[area_size, area_size]], [[int(45 - area_size/2), int(44 - area_size/2)]]
    step = 2

    # Place targets
    target_mesh_points = []
    for i in range(len(area_origin)):
        target_mesh_points = make_target_mesh(target_mesh_points,terrain,area_origin[i], area_dim[i], step)

    # Create target list
    target_localized_successfully = [1 for _ in range(len(target_mesh_points ))]
    targets_in = []
    for target in target_mesh_points:
        for pt in target:
            targets_in.append(pt) 

    # Define the list of optimizers being tested
    optimizers = ['DIRECT']

    ### RUN THE OPTIMIZER
    for i in range(2): # run 2x - first is trace, second is det!
        # Generate figs and plots
        fig = plt.figure()
        ax_opt = fig.add_subplot()
        step_num = i+1

        # Loop over the tested optimizers
        for optimizer in optimizers: 

            # For now, always consider a zero obj constraint
            obj_constraint = 0 
            
            # Comment out this section to optimize w/o the constraint
            #if i == 1: #second pass
                # Use constraint from ASSOCIATED optimizer!
            #    if optimizer == 'DIRECT':
            #        obj_constraint = obj_constraint_DIRECT
            #    elif optimizer == 'DE':
            #        obj_constraint = obj_constraint_DE
            #    elif optimizer == 'SA':
            #        obj_constraint = obj_constraint_SA

            # Construct tuples to pass in
            sensor_chars = (sensor_rad, sensor_type, num_sensors, sensor_comm_ratio, meas_type, obj_constraint, step_num)
            target_ins = (target_mesh_points, terrain)

            # Set the bounds (controls the number of variables to reason over)
            bounds = [(0, terrain_width-1),(0, terrain_height-1)]*(num_sensors) #bounds placing points to only on the map
            sensor_list = (sensor_chars, target_ins)   

            fcn_eval_list, fcn_counter = [], []
            counter, best_fcn = 0, 0
            
            print("Start", optimizer, "optimizer at iteration",i+1,"under area",area_size)
            print("counter, best_fcn, tr_sum, det_mult, sens_sens_penalty, valid_WSN_penalty, sens_target_penalty, valid_place_penalty, eps_opt_penalty")
            
            ### OPTIMIZER TYPES
            if optimizer == 'DIRECT':
                res = optimize.direct(objective_fcn, bounds=bounds, args = sensor_list, vol_tol=vol_tol, maxfun = maxfun_1, maxiter = max_iters)
                data_name_list.append("DIRECT Optimizer curve")
                print("Complete optimizing under", optimizer, "with value and iterations", res.fun, counter)
                print(res.message)
                
                x_out = res.x #save off prev results to start fine-tuning
                res = optimize.minimize(objective_fcn, x_out, args = sensor_list, method='BFGS', jac='3-point', options={'gtol': 0.01})
                obj_constraint_DIRECT = -1.0*res.fun
                print("Complete fine optimizing under", optimizer, "with value and iterations", res.fun, counter)

            elif optimizer == 'DE':
                res = optimize.differential_evolution(objective_fcn, bounds=bounds, args=sensor_list, strategy='best1bin', 
                                                    popsize=popsize, tol=0.01, mutation=mutation, recombination=recombination, maxiter=(int(max_fun_DE/(2*popsize*num_sensors))-1))
                data_name_list.append("Diff Ev Optimizer curve")
                print("Complete optimizing under", optimizer, "with value and iterations", res.fun, counter)
                print(res.message)
                
                x_out = res.x #save off prev results to start fine-tuning
                res = optimize.minimize(objective_fcn, x_out, args = sensor_list, method='BFGS', jac='3-point', options={'gtol': 0.01})
                obj_constraint_DE = -1.0*res.fun
                print("Complete fine optimizing under", optimizer, "with value and iterations", res.fun, counter)

            elif optimizer == 'SA':
                res = optimize.dual_annealing(objective_fcn, bounds=bounds, args=sensor_list, maxfun=maxfun_1, maxiter=max_iter_SA, no_local_search=True, 
                                            initial_temp=initial_temp, restart_temp_ratio=restart_temp_ratio, visit=visit, accept = -0.01)
                data_name_list.append("Sim Anneal Optimizer curve")
                print("Complete optimizing under", optimizer, "with value and iterations", res.fun, counter)
                print(res.message)
                
                x_out = res.x #save off prev results to start fine-tuning
                res = optimize.minimize(objective_fcn, x_out, args = sensor_list, method='BFGS', jac='3-point', options={'gtol': 0.01})
                obj_constraint_SA = -1.0*res.fun
                print("Complete fine optimizing under", optimizer, "with value and iterations", res.fun, counter)

            ### RESULTS PER RUN
            ax_opt.plot(fcn_counter, fcn_eval_list, linestyle='None', marker='+', label = optimizer)
            x_out = res.x

            # Save off results 
            sensor_positions_save.append(x_out)
            sensor_names_save.append("length"+str(area_size))

            inputs = target_localized_successfully, targets_in, x_out, sensor_rad, meas_type, terrain, 0 #LOS_flag == 1
            (FIMs_run, det_sum) = build_map_FIMS(inputs)
            det_mult1, tr_sum1 = return_obj_vals(FIMs_run)

            tr_vals.append(tr_sum1)
            det_vals.append(det_sum)
            FIM_areas, maj_axis, std_devs = [], [], []

            for q in range(len(targets_in)//2):
                target_i = [targets_in[0+2*q], targets_in[1+2*q]]

                # Compute the stats
                Cov_Mat, a, b, sx, sy, area_first = uncertainty_ellipse_stats(FIMs_run[q], target, 2.48, 1)
                FIM_areas.append(area_first)
                maj_axis.append(a)
                std_devs.append(sx**2)
                std_devs.append(sy**2)
            
            print("Mean area, largest axis, iter #:", np.mean(FIM_areas), np.max(maj_axis), i)
            print("Sensor positions:", x_out)

            if i == 0:
                tr_area_list.append(np.mean(FIM_areas))
                tr_largest_axis.append(np.max(maj_axis))
                tr_mean_var.append(np.mean(std_devs))
                opt_var.append("Trace")

            if i == 1:
                det_area_list.append(np.mean(FIM_areas))
                det_largest_axis.append(np.max(maj_axis))
                det_mean_var.append(np.mean(std_devs))
                opt_var.append("Det")

        # End for optimizer ... loop

    print("------------------------") #seperate out different areas

# Print the values and inspect
print("Trace solutions (first opt), then det solts (second)")
print(det_area_list, tr_area_list)
print(det_largest_axis, tr_largest_axis)
print(det_mean_var, tr_mean_var)
print(det_vals, tr_vals)
print(opt_var)

# Save the data off
area_save_off = [tr_area_list, det_area_list]
area_save_names = ["Trace optimization mean area","Det optimization mean area"]
save_data(area_save_off, area_save_names, "Area_vary_optimization_method", "Optimal solutions for a varied area under trace and det - area")

axis_save_off = [tr_largest_axis, det_largest_axis]
axis_save_names = ["Trace optimization max axis","Det optimization max axis"]
save_data(axis_save_off, axis_save_names , "axis_vary_optimization_method", "Optimal solutions for a varied area under trace and det - max axis")

var_save_off = [tr_mean_var, det_mean_var]
var_save_names = ["Trace optimization mean variance","Det optimization mean variance"]
save_data(var_save_off, var_save_names, "variance_vary_optimization_method", "Optimal solutions for a varied area under trace and det - mean variance")

var_save_off = [tr_vals, det_vals]
var_save_names = ["Trace vals","Det vals"]
save_data(var_save_off, var_save_names, "area_variance_tests_obj_vals", "Objective values for trace and det optimization / area test")

save_data([area_list], ["Tested area lengths"], "tested_area_lengths", "tested areas by length of side edge")

save_data(sensor_positions_save, sensor_names_save, "sensor_pos_save", "outputted sensor positions")


fig = plt.figure()
ax_areas = fig.add_subplot()
ax_areas.plot(area_plot, det_area_list, label='det area list')
ax_areas.plot(area_plot, tr_area_list, label='tr area list')
plt.legend()

plt.show()

# ---------PLOTTING AND VISUALIZATION------------------
# Comment out the random ones if we don't want to visualize them! 
# It's not always necessary because the ellipses can get unecessarily large

sensor_save_data, sensor_save_names = [], []

_map = make_basic_seismic_map(num_sensors, sensor_rad, sensor_type, meas_type, sensor_comm_ratio, x_out)
ax = terrain.plot_grid(_map)
new_map = add_targets(ax, targets_in)


plt.show()
