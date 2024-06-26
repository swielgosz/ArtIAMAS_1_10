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
def min_distance_penalty(target, sensor, d_min):
    '''
    ACCEPTS: target (2x1 list of floats), sensor (2x1 list floats), d_min (scalar of minimum distnace)
    RETURNS: penalty scalar
    '''
    mu = 5. # parameter in the sigmoid fcn
    tx, ty = target[0], target[1]
    
    # Loop over each sensor and add the penalty

    # distance computation
    sx, sy = sensor[0], sensor[1]
    d_ts = ((tx-sx)**2+(ty-sy)**2)**(1/2)
    
    # Return penalty
    return 1/(1+np.exp(-mu*(d_ts-d_min)))

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
    #return 1/(1+np.exp(-mu*(dist-d_min)))
    return 1

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

    # Return the average
    return (penalty_mult/(10*len(sensors)//2))
