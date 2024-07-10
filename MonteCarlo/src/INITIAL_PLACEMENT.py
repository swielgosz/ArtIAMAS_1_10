from helper_calculations.landscape import Configuration, Node
from helper_calculations.sensors import Sensor
from helper_calculations.utils import get_csv_size
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

# -------- DESCRIPTION -----------------
# This script does one iteration of the pareto front plot by finding the 
# optimum between localization and coverage using a tunable weighted value
# --------------------------------------

# --------------PARAMETERS------------------
# SCRIPT CONTROL
MAX_ITERS = 100000 # max number of iterations for Monte Carlo simulation
optimization_choice = 2 #1 = monte carlo, 2 = sim annealing
CASE = 1
SUBCASE = 2 #1 = sequential, 2 = complete (ex, optimize all sensors at once)
w = 0.999 # Weighting towards coverage vs localization
# TERRAIN INPUTS

my_path = Path(__file__).parent / "../data/terrain/GIS_terrain_100_100.csv"
rows, cols = get_csv_size(my_path)
terrain_width = cols
terrain_height = rows
terrain = Configuration(cols, rows)
terrain.load_from_csv(my_path)

# OPTIMIZER PARAMETERS
# https://docs.scipy.org/doc/scipy/reference/optimize.html
best_score = 0
counter = 0
best_config = None
config_history = []

# Sim annealing parameters
restart_temp_ratio, visit, no_local_search = 0.01, 1.25, True

# DIRECT parameters
#vol_tol = [1e-4, 1e-6, 1e-8, 1e-10, 1e-12, 1e-14, 1e-16]
vol_tol = [1e-24]

# SENSOR LIST
sensor_rad = [15, 15, 15, 20, 20]
sensor_type = ["acoustic","acoustic", "acoustic","seismic","seismic"]
num_sensors = len(sensor_type)
sensor_comm_ratio = 2 # ratio of sensor communication to sensing radius 
sensor_mode = ["radius","radius", "radius","bearing","bearing"] #radius or bearing
# ------------------------------------------

# THESE ARE SNIPPETS TRYING TO FIX ERROR RELATED TO MONTE CARLO GETTING STUCK IN LOOP DUE TO SPAMMING NETWORK CONNECTION
# #[ ...Snip... ]
# import smtplib
# #[ ...Snip... ]
# for USER in open(opts.USERS,'r'):
#     smtpserver = smtplib.SMTP(HOST,PORT)
#     smtpserver.ehlo()
#     verifyuser = smtpserver.verify(USER)
#     print("%s %s:  %s") % (HOST.rstrip(), USER.rstrip(), verifyuser)
#     smtpserver.quit()

# def online_check():
#   while True:
#     try:
#       con = urllib2.urlopen("http://www.google.com/")
#       data = con.read()
#       logging.debug('{0} Reached the host. Exiting online_check'.format(time.strftime('[ %H:%M:%S ]')))
#       break
#     except:
#       logging.debug('{0} Could not reach host trying again in 3 seconds'.format(time.strftime('[ %H:%M:%S ]')))
#       time.sleep(3)

def objective_fcn(x, *args):
    # unpack sensor params, initialize valid placements
    sensor_rad, sensor_type, num_sensors, sensor_comm_ratio, map_list = args
    valid_placement_check, valid_WSN = True, True
    
    global counter
    global fcn_eval_list
    global fcn_counter
    global SUBCASE

    counter += 1
    if counter % 50 == 0:
        print(counter, x)

    for check in range(len(x)//2):
        j, i = int(x[0+2*check]), int(x[1+2*check])

        # verify that the sensor is placed in a valid location
        if not terrain.is_valid_sensor_location(i, j):
            valid_placement_check = False
    
    # verify the sensor map is in a WSN
    # The LAST TWO POINTS ARE THE DESIGN POINTS IN SEQUENTIAL
    if SUBCASE == 1:
        map_list[-2], map_list[-1] = i, j 
        _map = make_basic_seismic_map(num_sensors, sensor_rad, sensor_type, sensor_mode, sensor_comm_ratio, map_list)
        valid_WSN = terrain.is_configuration_valid(_map)

    # Otherwise, all points are simutaneously optimized
    if SUBCASE == 2:
        map_list = x
        _map = make_basic_seismic_map(num_sensors, sensor_rad, sensor_type, sensor_mode, sensor_comm_ratio, map_list)
        valid_WSN = terrain.is_configuration_valid(_map)

    # Finally, perform the objective fcn computations
    if valid_placement_check and valid_WSN:
        # must return the NEGATIVE because we're minimizing
        coverage = terrain.get_configuration_observability(_map)
        localization = terrain.get_localization_score(_map)
        score = -(w*coverage + (1-w)*localization) #  get_configuration_observability will prioritize maximum coverage, get_localization_score will prioritize overlapping coverage
        fcn_eval_list.append(-1*score)
        fcn_counter.append(counter)
        if counter % 50 == 0:
            print(localization)
            print(coverage)

        return score
    
    # An invalid config will just return an arbitrarily large number
    else:
        score = 0
        fcn_eval_list.append(score)
        fcn_counter.append(counter)
        return score

try:  
    
    # Sequential DIRECT

    if optimization_choice == 2 and SUBCASE == 1:
        print('out')

        # sweep through the sensor arrangements
        for k in range(num_sensors):
            # set an initial condition
            fcn_eval_list, fcn_counter = [], []
            x0 = []
            i, j = np.random.randint(0, terrain_width), np.random.randint(0, terrain_height)
            x0.extend([i, j])

            # First sensor placed!
            if k ==  0:
                # inits for the optimizer
                sensor_placed_list = x0
                sensor_lists = (sensor_rad, sensor_type, num_sensors, sensor_comm_ratio, sensor_placed_list)
                bounds = [(0, terrain_width-1),(0, terrain_height-1)] #bounds placing points to only on the map
                
                # call the optimizer
                res = optimize.direct(objective_fcn, bounds=bounds, args = sensor_lists, vol_tol = vol_tol)
                #res = optimize.dual_annealing(objective_fcn, bounds=bounds, args = sensor_lists, 
                #                              restart_temp_ratio=restart_temp_ratio, visit = visit, no_local_search=no_local_search, maxfun=1000)
                
                # ouputs - overwrite points to list and plot objective fcn vs iteration 
                sensor_placed_list[-2], sensor_placed_list[-1] = res.x[0], res.x[1]
                print(res.message)
                fig = plt.figure()
                ax = fig.add_subplot()
                ax.plot(fcn_counter, fcn_eval_list, linestyle='None', marker='+')
            
            # Subsequent sensors placed
            else:
                # First compute the new limits - limits are the box that exactly encloses the sensor network
                x_sensor = [sensor_placed_list[2*i] for i in range(len(sensor_placed_list)//2)]
                y_sensor = [sensor_placed_list[2*i+1] for i in range(len(sensor_placed_list)//2)]

                # Use max, min to ensure that the limit does not exceed the map!
                xlower, xupper = max(min(x_sensor)-10*sensor_comm_ratio, 0), min(max(x_sensor)+10*sensor_comm_ratio, terrain_width)
                ylower, yupper = max(min(y_sensor)-10*sensor_comm_ratio, 0), min(max(y_sensor)+10*sensor_comm_ratio, terrain_height)

                # Now add the next sensor to optimize to the list 
                sensor_placed_list.extend([i, j])
                sensor_lists = (sensor_rad, sensor_type, num_sensors, sensor_comm_ratio, sensor_placed_list)
                bounds = [(xlower, xupper),(ylower, yupper)] #bounds placing points to only on the map
                
                # Optimize!
                res = optimize.direct(objective_fcn, bounds=bounds, args = sensor_lists, vol_tol = vol_tol)
                
                # Inpsect ouputs and plots
                print(res.message)
                sensor_placed_list[-2], sensor_placed_list[-1] = res.x[0], res.x[1]
                ax.plot(fcn_counter, fcn_eval_list,linestyle='None', marker='+')

            # Final Map
            _map = make_basic_seismic_map(num_sensors, sensor_rad, sensor_type, sensor_mode, sensor_comm_ratio, sensor_placed_list)   
            print('sensor', k+1, 'placed')
        
        x_final = sensor_placed_list

        # Customize some plots
        ax.set_ylabel('Objective Function Score'), ax.set_xlabel('Iteration Count')
        ax.set_ylim([-1, max(fcn_eval_list)*1.25])
        ax.set_title('Function Evaluation across Optimization')
        plt.show()

    # Complete DIRECT
    if optimization_choice == 2 and SUBCASE == 2:

        fig = plt.figure()
        ax = fig.add_subplot()
        x_res_plot_list = []

        for count, vol_tol_test in enumerate(vol_tol):
        
            for k in range(num_sensors):
                # set an initial condition
                fcn_eval_list, fcn_counter = [], []
                x0 = []
                i, j = np.random.randint(0, terrain_width), np.random.randint(0, terrain_height)
                x0.extend([i, j])

            # inits for the optimizer
            sensor_placed_list = x0
            sensor_lists = (sensor_rad, sensor_type, num_sensors, sensor_comm_ratio, sensor_placed_list)
            bounds = [(0, terrain_width-1),(0, terrain_height-1)]*(num_sensors) #bounds placing points to only on the map
                    
            # call the optimizer
            res = optimize.direct(objective_fcn, bounds=bounds, args = sensor_lists, vol_tol = vol_tol_test)

            print(res.message)
            ax.plot(fcn_counter, fcn_eval_list, linestyle='None', marker='+')
            ax.set_ylabel('Objective Function Evaluations'), ax.set_xlabel('Iteration Count')
            ax.set_ylim([-1, max(fcn_eval_list)*1.25])
            ax.set_title('Function Evaluation across Optimization')
            ax.grid(True, which='minor')   

            x_final = res.x
            x_res_plot_list.append(x_final)
            print('\nFINAL LOCATIONS:', x_final, 'TOLERANCE: ', vol_tol_test)

        fig = plt.figure()
        ax = fig.add_subplot()

        for i, list_plot in enumerate(x_res_plot_list):
            x_sensor = [list_plot[2*i] for i in range(len(list_plot)//2)]
            y_sensor = [list_plot[2*i+1] for i in range(len(list_plot)//2)]    
            ax.scatter(x_sensor, y_sensor)

        ax.set_ylabel('Map y'), ax.set_xlabel('Map x')
        ax.set_ylim([30, 70]), ax.set_xlim([30, 70])
        ax.set_title('Map of placed sensors')   
        ax.grid(True, which='minor')  
        plt.show()

    # MONTE CARLO OPTIMIZATION
    already_placed = []
    if optimization_choice == 1:
        # Run thru max iterations
        for i in range(MAX_ITERS):
            print("Iteration: ", i, end="\r")
            x = []
            for k in range(num_sensors):
                j = np.random.randint(0, terrain_width)
                i = np.random.randint(0, terrain_height)
            
                while (i, j) in already_placed or not terrain.is_valid_sensor_location(i, j):
                    j = np.random.randint(0, terrain_width)
                    i = np.random.randint(0, terrain_height)

                x.append(i)
                x.append(j)

            _map = make_basic_seismic_map(num_sensors, sensor_rad, sensor_type, sensor_comm_ratio, x)
            
            already_placed.append((i, j))

            # Calculate cost and check if configuration is valid
            score = terrain.get_configaration_observability(_map)
            valid = terrain.is_configuration_valid(_map)

            # Record each valid configuration that has a cost greater than the previous max
            if score > best_score and valid:
                print("New best score: ", score)
                best_score = score
                best_config = _map
                config_history.append(_map)
                x_final = x


    # Create a final plot
    _map = make_basic_seismic_map(num_sensors, sensor_rad, sensor_type, sensor_mode, sensor_comm_ratio, x_final)
    score = terrain.get_configuration_observability(_map)
    best_score = score
    best_config = _map
    ax = terrain.plot_grid(_map)
    print(best_score)
    pp.pprint(best_config)
    plt.show()
            
    # Save each new best configuration
    dirpath = './case_' + str(CASE)
    if os.path.exists(dirpath) and os.path.isdir(dirpath):
        shutil.rmtree(dirpath)
        os.mkdir(dirpath)
    else:
        os.mkdir(dirpath)
    terrain.save_history(config_history,dirpath)
    terrain.plot_grid(best_config)

# Save current best configuration if user stops before max iteration
except KeyboardInterrupt:
    print("Stopped early")
    print(best_score)
    pp.pprint(best_config)
    dirpath = './case_' + str(CASE) 
    if os.path.exists(dirpath) and os.path.isdir(dirpath):
        shutil.rmtree(dirpath)
        os.mkdir(dirpath)
    else:
        os.mkdir(dirpath)
    terrain.save_history(config_history,dirpath)
    terrain.plot_grid(best_config)