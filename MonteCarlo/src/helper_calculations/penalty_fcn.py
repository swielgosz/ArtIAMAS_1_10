import numpy as np
from .landscape import Configuration, Node

# This function code should accept some map and return the penalty associated with a constraint violation

# List of constraints we need to implement (? denotes a maybe):
    # WSN enforcement (equality)
    
    
    
    # Anything else???

# Equality constraint penalties are usually framed as a P = f(x)
# Inequality constraints are usually framed as P = if{x>limit}{0, g(x)} 
    # where it's zero if a constraint is satisfied and g(x) otherwise
    # DANIELA DE PALMA et al (2024) use a sigmoid to represent this, but we obviously dont have to

# Minimum sensor-target distance (inequality)
def min_distance_penalty(targets, sensors, d_min):
    '''
    ACCEPTS: target (2x1 list of floats), sensor (2x1 list floats), d_min (scalar of minimum distnace)
    RETURNS: penalty scalar
    '''
    sens_target_min = 100
    mu = 4. # parameter in the sigmoid fcn
    for i in range(len(targets)//2):
        for j in range(len(sensors)//2):
            tx, ty = targets[0+2*i], targets[1+2*i]
            sx, sy = sensors[0+2*j], sensors[1+2*j]
            d_ts = ((tx-sx)**2+(ty-sy)**2)**(1/2)
            if d_ts < sens_target_min:
                sens_target_min = d_ts
    
    # Return penalty
    return 1/(1+np.exp(mu*(sens_target_min-d_min)))#*((sens_target_min-d_min) + 1)
    #if sens_target_min < d_min:
    #    return 0.25*(sens_target_min - d_min)**2
    #else: 
    #    return 0


# Minimum sensor-sensor spacing (inequality)
def min_sensor_distance_penalty(sensors, d_min):
    '''
    ACCEPTS: target sensors (2Nx1 list floats, every 2 pairs corresponds to list of sensor pos)
    RETURNS: penalty scalar
    '''
    mu = 3. # parameter in the sigmoid fcn
    min_dist = 110.

    for i in range(len(sensors)//2): #check i'th sensor against all other
        for j in range(i+1, len(sensors)//2): # only check sensor pairs not yet checked
            sens1_x, sens1_y = sensors[0+2*i], sensors[1+2*i]
            sens2_x, sens2_y = sensors[0+2*j], sensors[1+2*j]
            dist = ((sens1_x-sens2_x)**2+(sens1_y-sens2_y)**2)**(1/2)

            # Check min condition
            if dist < min_dist:
                min_dist = dist

    # computes the penalty from the minimum distance of all sensor-sensor pairs
    return 1/(1+np.exp(-mu*(min_dist-d_min)))
    #return 1

# Valid placement (equality)
def valid_sensor_penalty(sensors, terrain):
# I'm doing a kinda kernel smoothing type thing: 
# https://en.wikipedia.org/wiki/Kernel_(image_processing)

    penalty_mult = 0.
    for i in range(len(sensors)//2): 
        sx, sy = round(sensors[0+2*i]), round(sensors[1+2*i])
        # Check main square
        if terrain.is_valid_sensor_location(sy, sx):
            penalty_mult += 2
        else:
            penalty_mult += 0
        
        # Check if sensor is at/near the edge of the plot
        # If so, return just 0 (worst-case penalty)
        if (sx < 98 and sy < 98 and sx > 1 and sy > 1):
            # Check all the plus/minus
            if terrain.is_valid_sensor_location(sy+1, sx):
                penalty_mult += 1
            else:
                penalty_mult += 0
            
            if terrain.is_valid_sensor_location(sy-1, sx):
                penalty_mult += 1
            else:
                penalty_mult += 0

            if terrain.is_valid_sensor_location(sy+2, sx):
                penalty_mult += 1
            else:
                penalty_mult += 0
            
            if terrain.is_valid_sensor_location(sy-2, sx):
                penalty_mult += 1
            else:
                penalty_mult += 0
            
            if terrain.is_valid_sensor_location(sy, sx+1):
                penalty_mult += 1
            else:
                penalty_mult += 0
            
            if terrain.is_valid_sensor_location(sy, sx-1):
                penalty_mult += 1
            else:
                penalty_mult += 0

            if terrain.is_valid_sensor_location(sy, sx+2):
                penalty_mult += 1
            else:
                penalty_mult += 0
            
            if terrain.is_valid_sensor_location(sy, sx-2):
                penalty_mult += 1
            else:
                penalty_mult += 0

        else:
            penalty_mult += 0

    # Return the average
    return (penalty_mult/(10*len(sensors)//2))

def WSN_penalty(sensors, comm_radii):
    sensors_list = []
    for i in range(len(sensors)//2):
        sensors_list.append([sensors[0+2*i], sensors[1+2*i]])

    # Just something large
    min_dist = 1000
    mu = 8

    # find the smallest value that violates the constraint
    for i, sensor_i in enumerate(sensors_list):
        for j, sensor_j in enumerate(sensors_list):
            if i == j:
                continue
            else:
                xi, yi = sensor_i[0], sensor_i[1]
                xj, yj = sensor_j[0], sensor_j[1]
                dist = ((xi-xj)**2+(yi-yj)**2)**(1/2)

                if ((dist <= min_dist) & (dist > comm_radii)):
                    min_dist = dist


    return 1/(1+np.exp(mu*(min_dist-comm_radii)))
