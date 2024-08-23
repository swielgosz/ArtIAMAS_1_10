import numpy as np
from .landscape import Configuration, Node
import random

def make_target_mesh(target_list, terrain, area_origin, area_dim, step):
    
    x_origin, y_origin = area_origin
    width, height = area_dim
   
    # Loop through rectangle and place targets
    for i in range(x_origin, x_origin + width, step): #changed these to step every 2 for TESTING - CHANGE BACK IF NECESSARY 
         for j in range(y_origin, y_origin + height, step):
              num_targets = terrain.get_target_density(j,i)
              terrain.grid[i][j].has_target = True
              terrain.grid[i][j].num_targets = num_targets # define the number of targets at each gridpoint
              #print(f"Number of targets at ({i}, {j}): {terrain.grid[i][j].num_targets}")
              for k in range(num_targets):
                   target_list.append([i,j]) # update list of targets

    return target_list               

def sample_target_mesh(mesh, n):
     return random.sample(mesh, n)

def upper_achieveable_limit(terrain, area_origins, area_dims, num_sens):
     # Loop through rectangle and place targets
     tot_targets = 0.

     # Calculate the total number of targets from mesh
     for kk in range(len(area_origins)): 
          area_origin, area_dim = area_origins[kk], area_dims[kk]

          x_origin, y_origin = area_origin
          width, height = area_dim

          for i in range(x_origin, x_origin + width, 2): #changed these to step every 2 for TESTING - CHANGE BACK IF NECESSARY 
               for j in range(y_origin, y_origin + height, 2):
                    num_targets = terrain.get_target_density(j,i)
                    terrain.grid[i][j].has_target = True
                    terrain.grid[i][j].num_targets = num_targets # define the number of targets at each gridpoint
                    # Compute the number of targets in mesh
                    tot_targets += num_targets
          
     print("Total number of targets: ", tot_targets)
     
     # when eta == 0, S = bearings only, and sigma = 1
     # trace and det upper limits
     return  (tot_targets*num_sens, tot_targets*num_sens**2/4)
