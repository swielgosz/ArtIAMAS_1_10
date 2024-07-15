import numpy as np
from .landscape import Configuration, Node


def make_target_mesh(target_list,terrain, area_origin, area_dim):
    
    x_origin, y_origin = area_origin
    width, height = area_dim
   
    # Loop through rectangle and place targets
    for i in range(x_origin, x_origin + width):
         for j in range(y_origin, y_origin + height):
              num_targets = terrain.get_target_density(i,j)
              terrain.grid[i][j].has_target = True
              terrain.grid[i][j].num_targets = num_targets # define the number of targets at each gridpoint
              print(f"Number of targets at ({i}, {j}): {terrain.grid[i][j].num_targets}")
              for k in range(num_targets):
                   target_list.append([i,j]) # update list of targets

    return target_list               
