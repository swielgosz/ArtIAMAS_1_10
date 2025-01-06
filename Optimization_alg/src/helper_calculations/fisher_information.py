import numpy as np
from .sensor_vision import get_LOS_coeff
from .landscape import Configuration, Node
from .penalty_fcn import min_distance_penalty

def build_FIM(sensors, target, sensor_types, sigma_e):
    '''
    Calculate the fisher information matrix for:
    1 point and
    N homogenous sensors that see a target
    about the localized target, p
    
    Returns the FIM following sources [1], [3]
    PENALTY CONSTRAINTS HARD-CODED HERE!!!
    '''

    # Penalty constraints
    d_min = 8

    # Need to first generate a list of phi's and r'i
    phi = []
    r = []

    # Construct a phi for each sensor relative to a single pt
    for i in range(len(sensors)//2):
        xi, yi = sensors[0+2*i], sensors[1+2*i]
        dx, dy = (xi- target[0]), (yi - target[1])
        ri = (dx**2+dy**2)**(1/2)
        if np.isclose(dx, 0, rtol=1e-06) == True:
            phi_i = np.pi
        else:
            phi_i = np.arctan2(dy, dx)

        phi.append(phi_i)
        r.append(ri)

    #print("phi list:", phi)
    #print("r list:", r)

    # Consturct the matrix
    FIM = np.array(([0.,0.],[0.,0.]))
    for i in range(len(phi)): #phi is one-to-one with sensor. Loops over sensor list effectively

        # Get penalty scalars
        # penalty = min_distance_penalty(target, [sensors[0+2*i], sensors[1+2*i]], d_min)
        penalty = 1 # NOTE: CHANGED TO GLOBALLY SATISFY CONSTRAINT

        # Distance-based noise scaling
        eta = 0.04
        dist_noise_scaling = (1+eta*r[i])

        if sensor_types[i] == "bearing":
            # If the sensor is placed on the target, avoid div by 0 error
            if np.isclose(r[i], 0, rtol=1e-06) == True:
                FIM += penalty*(1/(1)**2)*np.array(([(np.cos(phi[i]))**2,-0.5*np.sin(2*phi[i])],
                                [-0.5*np.sin(2*phi[i]),(np.sin(phi[i]))**2]))*(1/dist_noise_scaling)**2
            
            else: #note: should be 1/r^2, but I'm testing stuff - JM
                FIM += penalty*(1/(r[i])**2)*np.array(([(np.cos(phi[i]))**2,-0.5*np.sin(2*phi[i])],
                                [-0.5*np.sin(2*phi[i]),(np.sin(phi[i]))**2]))*(1/dist_noise_scaling)**2
            
            # Add in contributions from dist-dep noise
            FIM += 2*eta**2*penalty*np.array(([(np.sin(phi[i]))**2,0.5*np.sin(2*phi[i])],
                            [0.5*np.sin(2*phi[i]),(np.cos(phi[i]))**2]))*(1/dist_noise_scaling)**2

            # Add in the dist-noise scaling
            FIM = FIM

        if sensor_types[i] == "radial" or "radius":
            FIM += penalty*np.array(([(np.sin(phi[i]))**2,0.5*np.sin(2*phi[i])],
                            [0.5*np.sin(2*phi[i]),(np.cos(phi[i]))**2]))*((1+2*eta**2)/dist_noise_scaling)**2
        
            # Add in the dist-noise scaling
            # Add in contributions from dist-dep noise
            FIM = FIM
           
        FIM = FIM*(1/sigma_e[i])**2
        #print((1/sigma_e[i])**2, sigma_e[i])

    return FIM

def build_map_FIMS(inputs):
    '''
    This function constructs one FIM / target and returns a list of the numpy matrices
    ACCEPTS: target_localized_successfully - list of 0, 1's size N
             targets - list of target positions size 2N
             sensor_locs -  x,y sensor locations list size 2M
             sensor_rad - corresponding sensor radius size M
             meas_type - bearing or radial list
             terrain - the config node map
             LOS_flag - {0, 1} 1 if want to consider LOS, 0 otherwise
    RETURNS: list of numpy matrices, corresponding exactly to the number of targets
    '''

    target_localized_successfully, targets, sensor_locs, sensor_rad,meas_type, terrain, LOS_flag = inputs
    FIMs = []
    det_sum = 0.
    trace_sum = 0.
    for i in range(len(targets)//2):
        sub_sensor_FIM = []
        sub_sensor_type = []
        sub_sensor_sigma = []

        # For unlocalized targets
        if target_localized_successfully[i] == 0:
            FIM_p = np.array(([0, 0],[0,0]))
            FIMs.append(FIM_p)
            det_sum += np.linalg.det(FIMs[i])
            trace_sum += np.trace(FIMs[i])
        
        # For localized targets
        elif target_localized_successfully[i] == 1:
            xt, yt = targets[0+2*i], targets[1+2*i]
            
            # Check if the sensor sees the target
            for j in range(len(sensor_locs)//2):
                x, y = sensor_locs[0+2*j], sensor_locs[1+2*j] #Note, leave these as floats!
                dist = ((x-xt)**2+(y-yt)**2)**(1/2)
                if dist < sensor_rad[j]:
                    sub_sensor_FIM.append(x)
                    sub_sensor_FIM.append(y)
                    sub_sensor_type.append(meas_type[j])
                    if LOS_flag == 1:
                        sub_sensor_sigma.append(get_LOS_coeff([x,y], [xt, yt], terrain))
                    else:
                        sub_sensor_sigma.append(1)
            
            # Take the sub-list of sensors that see a target
            # and calculate the FIM
            # Penalty function is applied to EACH FIM!! - See build_FIM for application
            FIMs.append(build_FIM(sub_sensor_FIM, [xt, yt], sub_sensor_type, sub_sensor_sigma))
            det_sum += np.linalg.det(FIMs[i])
            trace_sum += np.trace(FIMs[i])
            
    # 2nd term is what we're returning to the objective!
    return FIMs, det_sum


def plot_uncertainty_ellipse(_map, FIM, target, confidence, plot_scale, color, solve_type):
    '''
    Accepts: map to plot, FIM cooresponding to ONE target, and one target
    Returns: plotted map w/ uncertainty ellipse plotted at target
    '''

    # Generate the covariance matrix
    Cov_Mat = np.linalg.inv(FIM)

    sx2 = (Cov_Mat[0,0])*1    #Sigma_x ^2 is the first element
    sy2 = (Cov_Mat[1,1])*1    #Sigma_y ^2 is the first element
    sxy = (Cov_Mat[0,1])*1    #Off-diagonal elements

    # Take sqrt to get sigma_x, sigma_y
    sx = sx2**(1/2)
    sy = sy2**(1/2)

    # Calculate a^2 and b^2, where a and b are the major axes
    a2 = (sx**2+sy**2)/2 + (((sx**2-sy**2)**2 )/ 4 + sxy**2)**(1/2)
    b2 = (sx**2+sy**2)/2 - (((sx**2-sy**2)**2 )/ 4 + sxy**2)**(1/2)

    # Get a and b, mult by the plot scale to correct for using different maps
    a = a2**(1/2)*plot_scale*1*confidence
    b = b2**(1/2)*plot_scale*1*confidence

    # Now calculate the rotation angle
    theta = (1/2)*np.arctan2(2*sxy, sx**2-sy**2)

    u=target[0]     #x-position of the center
    v=target[1]     #y-position of the center
    
    # Uncomment to print stats as necessary
    # print("ellipse area:", np.pi*a*b, "major axis:", a, "for target:", target)
    # print("Confid*sx, sy:", sx*plot_scale, sy*plot_scale, "major, minor axes", a, b, "for target:", target )

    # Plot!
    t = np.linspace(0, 2*np.pi, 100)
    if solve_type == 'numerical':
        linestyle = '--'
    elif solve_type == 'analytical':
        linestyle = '-'

    _map.plot(u+a*np.cos(t)*np.cos(theta) - b*np.sin(theta)*np.sin(t), 
              v+a*np.cos(t)*np.sin(theta) + b*np.cos(theta)*np.sin(t), 
              color = color, lw = 1.2, linestyle = linestyle)

    return _map


# fcn for computing the trace and det
def return_obj_vals(FIMs):
    det_mult = 1.
    tr_sum = 0.
    for kk in range(len(FIMs)):
        # FIM correspond to target list one-to-one
        det_mult = det_mult*np.linalg.det(FIMs[kk])
        tr_sum += np.trace(FIMs[kk])
    return det_mult, tr_sum

def uncertainty_ellipse_stats(FIM, target, confidence, plot_scale):
    '''
    Accepts: FIM cooresponding to ONE target, the target, confid level, plot scale
    Returns: stats associated with uncertainty ellipse as a tuple
             Returned tuple: (Cov matrix, a, b, sx, sy, area)
             a, b - major axes; sx, sy - std devs; area - uncertainty ellipse area at % confid
    '''

    # Generate the covariance matrix
    Cov_Mat = np.linalg.inv(FIM)

    sx2 = (Cov_Mat[0,0])*1    #Sigma_x ^2 is the first element
    sy2 = (Cov_Mat[1,1])*1    #Sigma_y ^2 is the first element
    sxy = (Cov_Mat[0,1])*1    #Off-diagonal elements

    # Take sqrt to get sigma_x, sigma_y
    sx = sx2**(1/2)
    sy = sy2**(1/2)

    # Calculate a^2 and b^2, where a and b are the major axes
    a2 = (sx**2+sy**2)/2 + (((sx**2-sy**2)**2 )/ 4 + sxy**2)**(1/2)
    b2 = (sx**2+sy**2)/2 - (((sx**2-sy**2)**2 )/ 4 + sxy**2)**(1/2)

    # Get a and b, mult by the plot scale to correct for using different maps
    a = a2**(1/2)*plot_scale*1*confidence
    b = b2**(1/2)*plot_scale*1*confidence

    # Finally, calculate the area
    area = np.pi*a*b

    return (Cov_Mat, a, b, sx, sy, area)