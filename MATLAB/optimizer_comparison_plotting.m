%% INTRODUCTION
% TITLE: Optimizer Evaluation Script
% PROJECT: Optimal heterogeneous WSN placements
% DATE: JUL 30, 2024
% AUTHORS: J. Mockler
% DESC: This script parses the optimizer outputted data and creates a plot

clear all; close all

%% READ IN DATA 
data_list = ["Case1_OptimizerRuns1.csv",
    "Case1_OptimizerRuns2.csv",
    "Case1_OptimizerRuns3.csv"
    "Case1_OptimizerRuns4.csv"];
n_samples = length(data_list);

% Initialize figures
figure (1)
hold on
figure (2)
hold on

SA_mean1 = zeros([100000, 1]);
DE_mean1 = zeros([100000, 1]);

SA_mean2 = zeros([100000, 1]);
DE_mean2 = zeros([100000, 1]);

% Plot all saved data runs!
for k = 1:length(data_list)

% Read in the list.
run = readtable(data_list(k), "NumHeaderLines", 0);
% If this is the first run, the DIRECT data is included, so we must use an
% extended member list
if k == 1
    mem_list1 = [1, 2, 3]; mem_list2 = [4, 5, 6];
else
    mem_list1 = [1, 2]; mem_list2 = [3, 4];
end

% Now loop through the csv file data
for j = 1:length(run{:,1})/2
    % pull name and data
    name_cell = run{2*(j-1)+1,1};
    name_str = name_cell{1,1};

    data_cell = run{2*(j-1)+1,2};
    data_str = data_cell{1,1};
    data = sscanf(data_str,'%f');
    
    % Pick appropriate fig to plot
    if ismember(j, mem_list1) % First set in each csv == first optimization
        figure (1)  
    elseif ismember(j, mem_list2) % Next set == second optimization
        figure (2)
    end
    
    % Pick the correct color and plot
    if strcmp(name_str, 'DIRECT Optimizer curve')
        % For DIRECT, we just need
        if ismember(j, mem_list1) % Plot DIRECT from first set
            figure (1)  
        elseif ismember(j, mem_list2)
            figure (2)
        end
        plot(1:length(data), data, 'color',"#0072BD", 'LineWidth', 1.5)
    
    elseif strcmp(name_str, 'Sim Anneal Optimizer curve')
        if ismember(j, mem_list1) % First three
            figure (1)
            % Add to list to get a mean profile at the end
            SA_mean1 = SA_mean1 + [data; zeros(100000-length(data), 1)];
        elseif ismember(j, mem_list2)
            figure (2)
            % Add to list to get a mean profile at the end
            SA_mean2 = SA_mean2 + [data; zeros(100000-length(data), 1)];
        end
        plot(1:length(data), data, 'color',"#A2142F", 'LineStyle','--')
    
    elseif strcmp(name_str, 'Diff Ev Optimizer curve')
        if ismember(j, mem_list1) % First three
            figure (1)
            DE_mean1 = DE_mean1 + [data; zeros(100000-length(data), 1)];
        elseif ismember(j, mem_list2)
            figure (2)
            DE_mean2 = DE_mean2 + [data; zeros(100000-length(data), 1)];
        end
        plot(1:length(data), data, 'color',"#EDB120", 'LineStyle','--')
    end

end % end csv parsing for statement

end % end data_list for statement

%% PLOT FORMATTING

figure (1) % first run fig
% Overlay the means
SA_mean1_cleaned = chop_array(SA_mean1);
DE_mean1_cleaned = chop_array(DE_mean1);
plot(1:length(SA_mean1_cleaned), SA_mean1_cleaned./n_samples, ...
    'color',"#A2142F", 'LineWidth', 1.5)
plot(1:length(DE_mean1_cleaned), DE_mean1_cleaned./n_samples, ...
    'color',"#EDB120", 'LineWidth', 1.5)

xlabel('Function Evaluation Count','Interpreter', 'latex') 
ylabel('Best Objective','Interpreter', 'latex')
title('Optimizer Runs - First Optimization','Interpreter', 'latex')
set(gca,'TickLabelInterpreter','latex')
legend('DIRECT', 'SA','DE','Interpreter',"latex")
grid minor
xlim([0, 12000])

figure (2) % sec run fig
SA_mean2_cleaned = chop_array(SA_mean2);
DE_mean2_cleaned = chop_array(DE_mean2);
plot(1:length(SA_mean2_cleaned), SA_mean2_cleaned./n_samples, ...
    'color',"#A2142F", 'LineWidth', 1.5)
plot(1:length(DE_mean2_cleaned), DE_mean2_cleaned./n_samples, ...
    'color',"#EDB120", 'LineWidth', 1.5)

xlabel('Function Evaluation Count','Interpreter', 'latex') 
ylabel('Best Objective','Interpreter', 'latex')
title('Optimizer Runs - Second Optimization','Interpreter', 'latex')
set(gca,'TickLabelInterpreter','latex')
legend('DIRECT', 'SA','DE','Interpreter',"latex")
grid minor
xlim([0, 12000])
