%% INTRODUCTION
% TITLE: Optimizer Evaluation Script
% PROJECT: Optimal heterogeneous WSN placements
% DATE: JUL 30, 2024
% AUTHORS: J. Mockler
% DESC: This script makes the computational load and analysis plot

clear; close all

%% Data
N = 2:1:8;

% Case 1 - data copy-and-pasted from simulation runs
DIRECT_wall_clock1 = [19.8	40.38	98.34	149.21]; % sec
SA_wall_clock1 = [216.79	321.448	547.479	617.67	799.4	789.95	1043.3];
DE_wall_clock1 = [66.3	221.41	797.92	1114.83	2373.314	3554.79 6313.38];

DIRECT_mean_area1 = [212.33	63.58	45.98	34.797];
SA_mean_area1 = [136.66	64.96	51.645	36.121	29.677	27.99	24.314];
DE_mean_area1 = [125.0166766	68.381	49.232	36.028	34.367	29.58	23.666];

DIRECT_mean_axis1 = [17.06	5.99	4.48	3.84];
SA_mean_axis1 = [10.246	6.151	5.171	4.054	3.459	3.492	3.361];
DE_mean_axis1 = [9.365	6.271	5.331	4.893	3.665	3.558	3.168];

% Case 2
DIRECT_wall_clock2 = [28.057	99.31	261.605	262.199	470.173 996.288	1145.746];
SA_wall_clock2 = [117.93 153.225	233.52	323.498	438.188	499.003 676.22];
DE_wall_clock2 = [23.061	110.46	260.978	630.243	961.9109	1340.796	1609.221];

DIRECT_mean_area2 = [107.46	45.802	37.77	26.399	24.952	19.747	17.025];
SA_mean_area2 = [101.659	56.432	46.862	29.213	23.37	23.903 17.256];
DE_mean_area2 = [135.873	54.614	37.496	27.219	26.668	18.559	16.874];

DIRECT_mean_axis2 = [9.851	4.878	4.46927371	3.558	3.521	2.984	2.707];
SA_mean_axis2 = [9.614	6.72	6.484	4.1414	3.759	4.132184284 3.074];
DE_mean_axis2 = [12.4	6.168	4.707	3.7	3.777516629	2.807	2.516];


%% Plotting

% Plots the wall-clock time
figure (1); hold on
plot(2:1:5, DIRECT_wall_clock1, 'o--', 'color',"#0072BD", 'LineWidth', 0.5 )
plot(N, SA_wall_clock1, 'o--','color',"#A2142F", 'LineWidth', 0.5)
plot(N, DE_wall_clock1 ,'o--','color',"#EDB120", 'LineWidth', 0.5)

plot(N, DIRECT_wall_clock2, '*--', 'color',"#0072BD", 'LineWidth', 0.5, 'MarkerFaceColor',"#0072BD", 'MarkerSize', 9 )
plot(N, SA_wall_clock2, '*--','color',"#A2142F", 'LineWidth', 0.5, 'MarkerFaceColor',"#A2142F",'MarkerSize',9 )
plot(N, DE_wall_clock2 ,'*--','color',"#EDB120", 'LineWidth', 0.5, 'MarkerFaceColor',"#EDB120", 'MarkerSize',9)

set(gca, 'YScale', 'log')
grid minor
ylabel('Optimization wall-clock time, sec', 'Interpreter','latex'); 
xlabel('$N$ sensors', 'Interpreter','latex')
legend('DIRECT, Case I', 'SA, Case I', ...
    'DE, Case I', 'DIRECT, Case II', 'SA, Case II', ...
    'DE, Case II','interpreter', 'latex')


% Plotting area
figure (2); 
subplot(2,1,1); hold on
plot(2:1:5, DIRECT_mean_area1, 'o--', 'color',"#0072BD", 'LineWidth', 0.5 )
plot(N, SA_mean_area1, 'o--','color',"#A2142F", 'LineWidth', 0.5)
plot(N, DE_mean_area1 ,'o--','color',"#EDB120", 'LineWidth', 0.5)

grid minor
ylabel('Mean area at 95\%', 'Interpreter','latex'); 
legend('DIRECT', 'SA', 'DE', 'interpreter', 'latex')
title('Case I: Uncertainty Ellipse Analysis', 'Interpreter','latex')

figure (3); 
subplot(2,1,1); hold on
plot(N, DIRECT_mean_area2, '*--', 'color',"#0072BD", 'LineWidth', 0.5, 'MarkerFaceColor',"#0072BD", 'MarkerSize', 9 )
plot(N, SA_mean_area2, '*--','color',"#A2142F", 'LineWidth', 0.5, 'MarkerFaceColor',"#A2142F",'MarkerSize',9 )
plot(N, DE_mean_area2 ,'*--','color',"#EDB120", 'LineWidth', 0.5, 'MarkerFaceColor',"#EDB120", 'MarkerSize',9)


%set(gca, 'YScale', 'log')
grid minor
ylabel('Mean area at 95\%', 'Interpreter','latex'); 
legend('DIRECT', 'SA', 'DE', 'interpreter', 'latex')
title('Case II: Uncertainty Ellipse Analysis', 'Interpreter','latex')

% Plotting mean major axis
figure (2)
subplot(2,1,2); hold on
plot(2:1:5, DIRECT_mean_axis1, 'o--', 'color',"#0072BD", 'LineWidth', 0.5 )
plot(N, SA_mean_axis1, 'o--','color',"#A2142F", 'LineWidth', 0.5)
plot(N, DE_mean_axis1 ,'o--','color',"#EDB120", 'LineWidth', 0.5)
grid minor
ylabel('Mean major axis at 95\%', 'Interpreter','latex')
xlabel('$N$ sensors', 'Interpreter','latex')

figure (3); 
subplot(2,1,2); hold on
plot(N, DIRECT_mean_axis2, '*--', 'color',"#0072BD", 'LineWidth', 0.5, 'MarkerFaceColor',"#0072BD", 'MarkerSize', 9 )
plot(N, SA_mean_axis2, '*--','color',"#A2142F", 'LineWidth', 0.5, 'MarkerFaceColor',"#A2142F",'MarkerSize',9 )
plot(N, DE_mean_axis2 ,'*--','color',"#EDB120", 'LineWidth', 0.5, 'MarkerFaceColor',"#EDB120", 'MarkerSize',9)

%set(gca, 'YScale', 'log')
grid minor
ylabel('Mean major axis at 95\%', 'Interpreter','latex'); 
xlabel('$N$ sensors', 'Interpreter','latex')
%legend('DIRECT, Case I', 'SA, Case I', ...
%     'DE, Case I', 'DIRECT, Case II', 'SA, Case II', ...
%     'DE, Case II','interpreter', 'latex')

%'color',"#0072BD", 'LineWidth', 1.5)
%'color',"#A2142F", 'LineStyle','--')
%'color',"#EDB120", 'LineStyle','--')