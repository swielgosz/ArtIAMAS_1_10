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
      
counter = 0
MAX_ITERS = 100000 # max number of iterations for Monte Carlo simulation
optimization_choice = 2
# Choose the optimization method - you can run both on this script
#1 = monte carlo, 2 = sim annealing
sensor_comm_ratio = 1.5 # ratio of sensor communication to sensing radius (because communication radius is typically larger)
CASE = 1

# Fcn accepts a list of sensor radii, corresponding type, and location
# The sensors in the list are placed on the map
def make_basic_seismic_map(num_sensors, sensor_rad, sensor_type, x):
    sensor_map = []
    already_placed = []

    for k in range(num_sensors):

        item = {}

        # pull data from lists
        item["j"] = x[1+2*k]
        item["i"] = x[0+2*k]
        item["sensor_type"] = sensor_type[k]
        item["sensor_radius"] = sensor_rad[k]
        item["sensor_comm_radius"] = sensor_comm_ratio*sensor_rad[k]

        # Record placed sensors
        sensor_map.append(item)
    
    # Return list of sensors
    return sensor_map


# Sim annealing objective function
# Fcn will return an arbitrarily large number if:
    #1. the sensor network is not fully connected
    #2. a sensor is placed on an invalid terrain
# Otherwise, the function will calculate and return the coverage score
# Per scipy.optimize requirements, it accepts a design vector x
# Here, x = [i_1, j_1, ... i_k, j_k] which are the coordinates of sensor k

def objective_fcn(x, *args):
    sensor_rad, sensor_type, num_sensors = args
    valid_placement_check = True
    # global counter
    # counter += 1
    # print(counter)
    
    # First check that a sensor is in a valid location
    for k in range(num_sensors):
        i = int(x[0+2*k])
        j = int(x[1+2*k])

        # if invalid, set to false 
        if not terrain.is_valid_sensor_location(i, j):
            valid_placement_check = False
    
    # Now check that we have a valid sensor network
    # Create a map of the sensor locations
    _map = make_basic_seismic_map(num_sensors, sensor_rad, sensor_type, x)
    valid_WSN = terrain.is_configuration_valid(_map)
    
    # Finally, perform the objective fcn computations
    if valid_placement_check and valid_WSN:
        # must return the NEGATIVE because we're minimizing
        global counter
        counter += 1
        print(counter)
        return -terrain.get_configuration_observability(_map)
        
    
    # An invalid config will just return an arbitrarily large number
    else:
        return 1000


# Define terrain 
# terrain_width = 1002
# terrain_height = 796
terrain_height = 100
terrain_width = 100
terrain = Configuration(terrain_width, terrain_height)
my_path = Path(__file__).parent / "../data/terrain/GIS_terrain_resize.csv"
terrain.load_from_csv(my_path)

# Initialize optimal sensor configuration cost and map
best_score = 0
best_config = None
config_history = []

# Define a list of sensor types and radii
# the design vector, x, will be a 2K length vector of (x, y) pairs for a sensor k
# the values are paired in order 
num_sensors = 4
sensor_rad = [10,10,10,10]
# sensor_rad = [5,5,5,5]
sensor_type = ["seismic", "acoustic", "seismic", "acoustic"]

try:  
    
    # SIM ANNEALING OPTIMIZATION
    if optimization_choice == 2:
        # First create a vector of initial condition
        x0 = []
        for k in range(num_sensors):
            j = np.random.randint(0, terrain_width)
            i = np.random.randint(0, terrain_height)
            x0.append(i)
            x0.append(j)

        # pass in the required lists to the objective
        sensor_lists = (sensor_rad, sensor_type, num_sensors)
        bounds = [(0, 99)]*2*num_sensors #bounds placing points to only on the map
        # Some documentation!!
        # https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.dual_annealing.html
        # res = optimize.dual_annealing(objective_fcn, bounds=bounds, args = sensor_lists)
        res = optimize.dual_annealing(objective_fcn, bounds=bounds, args = sensor_lists, no_local_search=True)
        # res = optimize.dual_annealing(objective_fcn, bounds=bounds, args = sensor_lists, restart_temp_ratio=0.01, visit = 1.25, no_local_search=True)
        # Convex minimization scipy fcn - may not actually work, I havent tried it in a while
        # res = minimize(objective_fcn, x0, args = sensor_lists)

        print('done')
    
        x_final = res.x

    # MONTE CARLO OPTIMIZATION
    already_placed = []
    if optimization_choice == 1:
        # Run thru max iterations
        for i in range(MAX_ITERS):
            time.sleep(0.01)
            if i%100 == 0:
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

            _map = make_basic_seismic_map(num_sensors, sensor_rad, sensor_type, x)
            
            already_placed.append((i, j))

            # Calculate cost and check if configuration is valid
            score = terrain.get_configuration_observability(_map)
            valid = terrain.is_configuration_valid(_map)

            # Record each valid configuration that has a cost greater than the previous max
            if score > best_score and valid:
                print("New best score: ", score)
                best_score = score
                best_config = _map
                config_history.append(_map)
                x_final = x


    # Create a final plot
    _map = make_basic_seismic_map(num_sensors, sensor_rad, sensor_type, x_final)
    score = terrain.get_configuration_observability(_map)
    best_score = score
    best_config = _map
    print(best_score)
    pp.pprint(best_config)
            
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