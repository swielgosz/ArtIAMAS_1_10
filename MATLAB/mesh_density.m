function [target_locs_list, density] = mesh_density(target_mesh)

    target_locs_list = [];
    density = [];
    density_counter = 1;

    for i = 1:length(target_mesh(:,1))
        if i == 1
            target_locs_list = [target_locs_list; target_mesh(1,:)];
        else
            if target_locs_list(end, :) == target_mesh(i, :)
                density_counter = density_counter + 1;
            else
                density = [density, density_counter];
                density_counter = 1;
                target_locs_list = [target_locs_list; target_mesh(i,:)];
            end
        end
    end

    density = [density, 1];

end