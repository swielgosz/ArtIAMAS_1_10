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

# -------- DESCRIPTION -----------------
#  Script is used to construct a pareto front between coverage vs
# localization in the initial placement. Of course, you can run just
# one point of the front (by picking the weighting fcn value)
# to get an optimal solution if just the sensor points are desired
# --------------------------------------

# --------------PARAMETERS------------------
# SCRIPT CONTROL
MAX_ITERS = 100000 # max number of iterations for Monte Carlo simulation
optimization_choice = 2 #1 = monte carlo, 2 = sim annealing
CASE = 1
SUBCASE = 2 #1 = sequential, 2 = complete (ex, optimize all sensors at once)
w = 0.999 # Weighting towards coverage vs localization
# TERRAIN INPUTS
terrain_height = 100 #796
terrain_width = 100 #1002
terrain = Configuration(terrain_width, terrain_height)
my_path = Path(__file__).parent / "../data/terrain/GIS_terrain_resize.csv"
terrain.load_from_csv(my_path)

# ----------OPTIMIZER PARAMETERS------------
# https://docs.scipy.org/doc/scipy/reference/optimize.html
best_score = 0
counter = 0
best_config = None
config_history = []

# DIRECT parameters
vol_tol = 1e-24

# SENSOR LIST
sensor_rad = [15, 15, 20, 20, 10]
sensor_type = ["acoustic","acoustic", "seismic","seismic","acoustic"]
num_sensors = len(sensor_type)
sensor_comm_ratio = 1.75 # ratio of sensor communication to sensing radius 
sensor_mode = ["radius","radius", "bearing", "bearing","radius"] #radius or bearing

# PARETO FRONT
w_vals = [0.155]
# ------------------------------------------


# ----------OBJECTIVE FUNCTION------------
def objective_fcn(x, *args):
    # unpack sensor params, initialize valid placements
    sensor_rad, sensor_type, num_sensors, sensor_comm_ratio, map_list, w = args
    valid_placement_check, valid_WSN = True, True
    
    global counter
    global fcn_eval_list
    global fcn_counter
    global SUBCASE

    counter += 1

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
        
        # Inspect progress
        #if counter % 200 == 0:
            #print('\nIter:',counter,'Localization score:',localization,'Coverage score:',coverage)
            #print('Positions:',x)

        return score
    
    # An invalid config will just return an arbitrarily large number
    else:
        score = 0
        fcn_eval_list.append(score)
        fcn_counter.append(counter)
        
        # Inspect progress
        #if counter % 200 == 0:
            #print('\nIter:',counter,'Localization score:',0,'Coverage score:',0)
            #print('Positions:',x)
        
        return score


# ---------- PARETO FRONT CONSTRUCTION --------------
coverage_evals, localization_evals = [], []
for kk in range(len(w_vals)):
    w_test = w_vals[kk]
    for k in range(num_sensors):
        # set an initial condition
        fcn_eval_list, fcn_counter = [], []
        x0 = []
        i, j = np.random.randint(0, terrain_width), np.random.randint(0, terrain_height)
        x0.extend([i, j])

    # inits for the optimizer
    sensor_placed_list = x0
    sensor_lists = (sensor_rad, sensor_type, num_sensors, sensor_comm_ratio, sensor_placed_list, w_test)
    bounds = [(0, terrain_width-1),(0, terrain_height-1)]*(num_sensors) #bounds placing points to only on the map

    # call the optimizer
    res = optimize.direct(objective_fcn, bounds=bounds, args = sensor_lists, vol_tol = vol_tol)

    print(res.message) 
    x_final = res.x
    
    # Get localization and coverage scores
    _map = make_basic_seismic_map(num_sensors, sensor_rad, sensor_type, sensor_mode, sensor_comm_ratio, x_final)
    coverage = terrain.get_configuration_observability(_map)
    localization = terrain.get_localization_score(_map)
    coverage_evals.append(coverage)
    localization_evals.append(localization)

    # Print and inspect
    print('\nFINAL LOCATIONS: ', x_final, 'TOLERANCE: ', vol_tol)
    print('\nWeighting:',w_test,'Coverage score:', coverage, 'Localization: ', localization)


# Create a final plot
print(coverage_evals)
print(localization_evals)
fig = plt.figure()
ax = fig.add_subplot()
ax.plot(localization_evals, coverage_evals, linestyle='None', marker='+')
ax.set_ylabel('Coverage Objective'), ax.set_xlabel('Localization Objective')
ax.set_title('Pareto Front')
plt.show()

# Create a final map for inspection
ax = terrain.plot_grid(_map)
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

