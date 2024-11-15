%% INTRODUCTION
% TITLE: Area variance testing plotting
% PROJECT: Optimal heterogeneous WSN placements
% DATE: 7 AUG 24
% AUTHORS: J. Mockler
% DESC: Plots the results of the area vs variance tests I ran. 

close all;

%% DATA
lengths = [24, 20, 16, 10, 4, 26, 22, 18, 14, 12, 6, 8];
% From EXPLORE optimization
avg_areas = [ 64.89087190992713 54.43954967885954 46.99278227294276 ...
    41.217183613994855 29.751320668615673 70.15441723156246 ...
    59.40521889843671 51.61587150287622 42.093433045795734 ...
    40.075776234861934 32.73348468181882 36.1801359508887];

% From EXPLOIT optimization
% represents the average of 4 trials, where 5 targets are randomly sampled
% from the associated square mesh
avg_areas_2nd = [34.73715774, 31.61915314, 30.13662651, 26.63287313, ...
    19.76826007, 38.67500323, 36.79789106, 31.2696338, 26.48511034, ...
    25.48984853, 20.78259159, 21.74922078];

Trace_vals=[116.42644315912122 91.92396040602847 67.63371342956825 27.0933725746942 5.540915143602352] ;
det_vals=[17.073801256958795 15.707888546280682 14.286095111941243 5.769609902379456 1.7153096279939561] ;
area_sub = [ 24 20 16 10 4 ];
largest_axis = [12.600268364254816 9.99730803561101 8.914230287466584 6.611224795348493 3.98036404221904];
avg_area_sub = [64.89087190992713 54.43954967885954 46.99278227294276 ...
    41.217183613994855 29.751320668615673];
trace_per_target = Trace_vals ./ (area_sub.^2);
det_per_target = det_vals ./ (area_sub.^2);

%% PLOTTING
figure (1)
scatter(lengths, avg_areas); hold on;
scatter(lengths, avg_areas_2nd)
grid minor
xlabel('Potential Mesh Area (side length of square)','Interpreter', 'latex') 
ylabel('Mean Uncertainty Ellipse Area at 95\% Confid.','Interpreter', 'latex')
set(gca,'TickLabelInterpreter','latex')
title('Effect of Initial Search Area on Performance','Interpreter', 'latex', 'Fontsize', 13)
legend('Explore Placement', 'Exploit Placement', 'location', 'northwest')

figure (2)
yyaxis left
scatter(avg_area_sub, trace_per_target); hold on;
ylabel('Average Trace per Target','Interpreter', 'latex') 
yyaxis right
scatter(avg_area_sub, det_per_target)
ylabel('Average Det per Target','Interpreter', 'latex') 

xlabel('Mean Optimal Uncert Ellipse Area','Interpreter', 'latex') 
set(gca,'TickLabelInterpreter','latex')
title('Explore Placement Optimal Values','Interpreter', 'latex', 'Fontsize', 13)
xlim([25, 70])
grid minor