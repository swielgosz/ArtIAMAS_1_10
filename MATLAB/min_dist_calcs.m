function [min_dist, min_target] = min_dist_calcs(sensor, target_locs_list)

min_dist = 100;
min_target = [0,0];
for i=1:length(target_locs_list(:,1))
    target = target_locs_list(i,:);
    dist = norm(target - sensor, 2);
    if dist < min_dist
        min_dist = dist;
        min_target = target;
    end
end

end