clear
close all
clc

%% Model Options
drug = "hydrophilic";
conditions =["human_avg","macaque_avg"];
isNonDimensionalized = true;

%phi_IT = 1;
%D_F = 5.7e-6;
%h_T = .3;
filePath = 'scaling_analysis.xlsx';
startRow = 1;
sheetArchetypal='ArchetypalComparisons';

M_0 = 30; %30mg LOADING

%cs_avg_conditions = cell(1,length(conditions));
%c0_conditions = cell(1,length(conditions));
%m0_conditions = cell(1,length(conditions));
T_conditions = cell(1,length(conditions));
t = linspace(0, 60*60*24*144, 10000);
t_day = t./(60*60*24);

%% ARCHETYPAL ANALYSIS

for i = 1:length(conditions)
    condition = conditions(i);
    [C, CI_avg_4C, CF_avg_4C, CE_avg_4C, CS_avg_4C, M_4C, N, Mass_remaining_4C, mR_inner, CB_4C, C_0, M_0, T, params] = solve_diffusion_5C(t, condition, drug, 'nonDimensional', true);

    startRow = writeExcelTable(T,filePath,sheetArchetypal,startRow,condition);
    %cs_avg_conditions{i} = CS_avg_4C;
    %c0_conditions{i} = C_0;
    %m0_conditions{i} = M_0;
    T_conditions{i} = T;

end

%% C*28/C0 calculation

C0_optimal_Cstar28 = T_conditions{1}.("C*28/C0")(4)/T_conditions{2}.("C*28/C0")(4);
result = sprintf("Optimal macaque C_0 given C*28/C0 = %.4f", C0_optimal_Cstar28);
resultRange = sprintf('A%d', startRow);
writecell({result}, filePath, 'Sheet', sheetArchetypal, 'Range', resultRange);

%% RMS Calculation

i28 = find(t >= 28*60^2*24,1,'first');
rms_MATLAB = rmse(cs_avg_conditions{1}(1:i28),cs_avg_conditions{2}(1:i28));
rms_CUSTOM = rmse_custom(t,cs_avg_conditions{1},cs_avg_conditions{2},'time_pt',28*24*3600);

%% RMS Minimization

drug = "hydrophilic";
condition = "human_avg";

% human calculation
[~, ~, ~, ~, cs_avg_human, ~, ~, ~, ~, ~, ~, ~, ~, ~] = solve_diffusion_5C(t, condition, drug, 'C_0', 1);

% Initial parameters for optimization
initialGuess_C0 = 1;
error_func_args = {'time_pt', 28*24*3600};
condition = "macaque_avg";

% Call the function with the error function arguments
[optimal_C0, C0_rmse_log, cs_avg_macaque_optimal] = runOptimization(t, cs_avg_human, drug, condition, @rmse_custom, initialGuess_C0, error_func_args);

%% Writing Optimization Results

sheetArchetypalOptimization='ArchetypalOptimizationMacaque';
startRow = 1;
C0_rmse_table = array2table(C0_rmse_log, 'VariableNames', {'C0_macaque', 'M0_macaque', 'RMSE'});
heading = sprintf("Optimal in Macaques: C_0 = %.2f", optimal_C0);

startRow = writeExcelTable(C0_rmse_table,filePath,sheetArchetypalOptimization,startRow,heading);


%% GETTING REMAINING DATA FOR PLOTTING COMPARISONS

drug = "hydrophilic";
condition = "macaque_avg";
[~, ~, ~, ~, allometric_optimal_cs_avg, ~, ~, ~, ~, ~, ~, ~, ~, ~] = solve_diffusion_5C(t, condition, drug, 'C_0', 1.04);
[~, ~, ~, ~, cstar28_optimal_cs_avg, ~, ~, ~, ~, ~, ~, ~, ~, ~] = solve_diffusion_5C(t, condition, drug, 'C_0', C0_optimal_Cstar28);
%[~, ~, ~, ~, rmse_optimal_cs_avg, ~, ~, ~, ~, ~, ~, ~, ~] = solve_diffusion_5C(t, condition, drug, 'C_0', optimal_C0);

rmse_optimal_cs_avg = cs_avg_macaque_optimal;

%% PLOTTING COMPARISONS

figure()
plot(t/(24*3600),cs_avg_human,'k--','LineWidth',4);
hold on
plot(t/(24*3600),cstar28_optimal_cs_avg,'LineWidth',4)
plot(t/(24*3600),rmse_optimal_cs_avg,'LineWidth',4)
plot(t/(24*3600),allometric_optimal_cs_avg,'LineWidth',4)

set(gca,'FontSize',28)
xlabel('Time (days)', 'FontSize', 36);
ylabel('Cs/C0','FontSize',36)
xlim([0,28])
legend({"Human", "Rule 1 - C*28/C0","Rule 2 - RMSE minimization", "Rule 3 - Allometric scaling"})


%% Parametric combinations - Running all conditions non-dimensionalized

drug = "hydrophilic";
isNonDimensionalized = true;
conditions =["human_avg","human_large","human_small","macaque_avg","macaque_large","macaque_small"];
filePath = 'scaling_analysis.xlsx';
startRow = 1;
sheetComparisonsC28Star='MultipleParametricComparisons';

human_conditions = conditions(contains(conditions, "human"));
%cs_avg_conditions = cell(1,length(conditions));
%c0_conditions = cell(1,length(conditions));
%m0_conditions = cell(1,length(conditions));
T_conditions = cell(1,length(conditions));
cs_avg_human_baselines = cell(1,length(human_conditions));

for i = 1:length(conditions)
    condition = conditions(i);
    [C, CI_avg_4C, CF_avg_4C, CE_avg_4C, CS_avg_4C, M_4C, N, Mass_remaining_4C, mR_inner, CB_4C, C_0, M_0, T, params] = solve_diffusion_5C(t, condition, drug, 'nonDimensional', true);
    
    startRow = writeExcelTable(T,filePath,sheetComparisonsC28Star,startRow,condition);
    %cs_avg_conditions{i} = CS_avg_4C;
    %c0_conditions{i} = C_0;
    %m0_conditions{i} = M_0;
    T_conditions{i} = T;

    % If the current condition is a human condition, store the cs_avg vector
    if contains(condition, "human")
        human_index = find(contains(human_conditions, condition));
        cs_avg_human_baselines{human_index} = CS_avg_4C;
    end

end

%% Parametric combinations - C28*/C0 calculation

% Extract human and macaque conditions directly
human_conditions = conditions(contains(conditions, "human"));
macaque_conditions = conditions(contains(conditions, "macaque"));

% Count the number of human and macaque conditions
numHumanConditions = length(human_conditions);
numMacaqueConditions = length(macaque_conditions);

% Identify human and macaque indices
human_indices = find(contains(conditions, "human"));
macaque_indices = find(contains(conditions, "macaque"));

% Calculate the total number of valid combinations
totalCombinations = numHumanConditions * numMacaqueConditions;

% Pre-allocate the comparisonResults cell array for efficiency
comparisonResults = cell(totalCombinations, 3); % 3 columns for 'Condition 1', 'Condition 2', and 'C0_optimal_Cstar28'
c0_optimal_cstar28 = cell(numHumanConditions, numMacaqueConditions);
cs_avg_cstar28 = cell(numHumanConditions, numMacaqueConditions);

% Initialize an index to fill in comparisonResults
resultIndex = 1;

% Iterate over each human and macaque condition
for i = 1:numHumanConditions
    for j = 1:numMacaqueConditions       
        % Get the indices for the current human and macaque conditions
        human_index = human_indices(i);
        macaque_index = macaque_indices(j);

        % Calculate the metric and run the model
        c0_optimal = T_conditions{human_index}.("C*28/C0")(4) / T_conditions{macaque_index}.("C*28/C0")(4);
        c0_optimal_cstar28{i,j} = c0_optimal;
        condition_macaque = conditions(macaque_index);
        [~, ~, ~, ~, cstar28_optimal_cs_avg, ~, ~, ~, ~, ~, ~, ~, ~, ~] = solve_diffusion_5C(t, condition_macaque, drug, 'C_0', c0_optimal);
        cs_avg_cstar28{i,j} = cstar28_optimal_cs_avg;

        % Store the result in the pre-allocated cell array
        comparisonResults{resultIndex, 1} = human_condition;
        comparisonResults{resultIndex, 2} = macaque_condition;
        comparisonResults{resultIndex, 3} = c0_optimal_cstar28{i,j};

        % Increment the index
        resultIndex = resultIndex + 1;
    end
end

% Convert the results to a table
comparisonTable = cell2table(comparisonResults, 'VariableNames', {'Condition1', 'Condition2', 'c0_optimal_cstar28'});

% Write the table to the Excel sheet
startRow = writeExcelTable(comparisonTable, filePath, sheetComparisonsC28Star, startRow, 'C*28/C0 Comparisons');

%% RMS Minimization for Multiple Combinations

drug = "hydrophilic";
%conditions = ["human_avg", "human_large", "human_small", "macaque_avg", "macaque_large", "macaque_small"];
filePath = 'scaling_analysis.xlsx';
startRow = 1;
sheetMultipleParameterOptimization = 'MultCombsOptimizationMacaque';

% Identify human and macaque indices
human_indices = find(contains(conditions, "human"));
macaque_indices = find(contains(conditions, "macaque"));

% Get the number of human and macaque conditions
numHumanConditions = length(human_indices);
numMacaqueConditions = length(macaque_indices);

% Initialize cell arrays with the correct size
cs_avg_optimal_conditions = cell(numHumanConditions, numMacaqueConditions);
C0_rmse_logs = cell(numHumanConditions, numMacaqueConditions);
optimal_C0_values = cell(numHumanConditions, numMacaqueConditions);

% Iterate over all human and macaque combinations
for i = 1:numHumanConditions
    for j = 1:numMacaqueConditions
        % Get the indices for the current human and macaque conditions
        human_index = human_indices(i);
        macaque_index = macaque_indices(j);

        % Perform human calculation for reference
        condition_human = conditions(human_index);
        [~, ~, ~, ~, cs_avg_human, ~, ~, ~, ~, ~, ~, ~, ~, ~] = solve_diffusion_5C(t, condition_human, drug, 'C_0', 1);

        % Set up optimization for macaque condition
        initialGuess_C0 = 1;
        error_func_args = {'time_pt', 28*24*3600};
        condition_macaque = conditions(macaque_index);

        % Run optimization for the current pair
        [optimal_C0, C0_rmse_log, cs_avg_macaque_optimal] = runOptimization(t, cs_avg_human, drug, condition_macaque, @rmse_custom, initialGuess_C0, error_func_args);

        % Store the results for plotting and analysis
        cs_avg_optimal_conditions{i, j} = cs_avg_macaque_optimal;
        C0_rmse_logs{i, j} = C0_rmse_log;
        optimal_C0_values{i, j} = [optimal_C0, C0_rmse_log(end, 2)]; % [C0_optimal, M0_optimal]

        % Create a table from the C0_rmse_log
        C0_rmse_table = array2table(C0_rmse_log, 'VariableNames', {'C0_macaque', 'M0_macaque', 'RMSE'});
        heading = sprintf("Optimal C0 for %s vs. %s: C0 = %.2f", condition_human, condition_macaque, optimal_C0);

        % Write the table to the Excel sheet
        startRow = writeExcelTable(C0_rmse_table, filePath, sheetMultipleParameterOptimization, startRow, heading);
    end
end

%% Writing Summary Results for RMSE optimization of multiple comparisons

% Preallocate the summaryResults cell array
summaryResults = cell(numHumanConditions * numMacaqueConditions, 4); % 4 columns for Condition1, Condition2, Optimal_C0, Optimal_M0

% Fill in the summaryResults with optimal values
resultIndex = 1;
for i = 1:numHumanConditions
    for j = 1:numMacaqueConditions
        % Store human and macaque conditions
        summaryResults{resultIndex, 1} = human_conditions{i};
        summaryResults{resultIndex, 2} = macaque_conditions{j};
        
        % Store the optimal C0 and M0 values
        summaryResults{resultIndex, 3} = optimal_C0_values{i, j}(1); % Optimal_C0
        summaryResults{resultIndex, 4} = optimal_C0_values{i, j}(2); % Optimal_M0
        
        % Increment the index
        resultIndex = resultIndex + 1;
    end
end

% Convert to table with appropriate column names
summaryTable = cell2table(summaryResults, 'VariableNames', {'Condition1', 'Condition2', 'Optimal_C0', 'Optimal_M0'});

% Write the summary table to the Excel sheet
startRow = writeExcelTable(summaryTable, filePath, sheetMultipleParameterOptimization, startRow, 'Summary of Optimal C0 and M0 for All Combinations');

%% RUNNING ALLOMETRIC ESTIMATIONS

% Define the allometric C0 values matrix (numHumanConditions x numMacaqueConditions)
allometric_C0_values = [
    1.04, 1.87, 0.63;  % human_avg vs macaque_avg, macaque_large, macaque_small
    0.65, 1.17, 0.39;  % human_large vs macaque_avg, macaque_large, macaque_small
    1.96, 3.52, 1.17   % human_small vs macaque_avg, macaque_large, macaque_small
];

cs_avg_allometric_conditions = cell(numHumanConditions, numMacaqueConditions);

for i = 1:length(human_indices)
        human_index = human_indices(i);
    for j = 1:length(macaque_indices)
        macaque_index = macaque_indices(j);

        C0_allometric = allometric_C0_values(i, j);
        [~, ~, ~, ~, cs_avg_allometric, ~, ~, ~, ~, ~, ~, ~, ~, ~] = solve_diffusion_5C(t, conditions(macaque_index), drug, 'C_0', C0_allometric);
        cs_avg_allometric_conditions{i,j} = cs_avg_allometric;
    end
end


%% SAVING .MAT FILE

% Define the conditions array
%conditions = ["human_avg", "human_large", "human_small", "macaque_avg", "macaque_large", "macaque_small"];

% Create separate arrays for human and macaque conditions
human_conditions = conditions(contains(conditions, "human"));
macaque_conditions = conditions(contains(conditions, "macaque"));

% Save all relevant variables for plotting in a .mat file
outputFileName = 'pk_metrics_data.mat';
save(outputFileName, 't', 'human_conditions', 'macaque_conditions', ...
    'cs_avg_human_baselines', ... % HUMAN BASELINES 
    'c0_optimal_cstar28','cs_avg_cstar28', ... %C%28 ANALYSIS
    'optimal_C0_values', 'cs_avg_optimal_conditions', ... %RMSE ANALYSIS
    'allometric_C0_values', 'cs_avg_allometric_conditions', ... %ALLOMETRIC ANALYSIS
    'drug', '-v7.3');


%% PLOTTING RESULTS

% Assuming you have already created the variables:
% t: time vector
% conditions: array of condition names (e.g., ["human_avg", "human_large", "human_small", "macaque_avg", "macaque_large", "macaque_small"])
% cs_avg_conditions: cell array containing Cs/C0 results for each condition
% cs_avg_optimal_conditions: cell array containing RMSE minimization results for each human-macaque pair
% allometric_C0_values: matrix of allometric C0 values (3x3 matrix for 3 human and 3 macaque conditions)
% solve_diffusion_5C: function handle to your model function
% drug: string representing the drug type (e.g., "hydrophilic")

%plotPKComparisons(t, conditions, cs_avg_human_baselines, cs_avg_cstar28, cs_avg_optimal_conditions, cs_avg_allometric_conditions);
plotPKVariability(t, conditions, cs_avg_human_baselines, cs_avg_cstar28, cs_avg_optimal_conditions, cs_avg_allometric_conditions)

%% LUTEAL ANALYSIS

%% Parametric combinations - Running all conditions non-dimensionalized

drug = "hydrophilic";
isNonDimensionalized = true;
conditions =["human_avg_luteal","human_large_luteal","macaque_avg_luteal","macaque_small_luteal"];
filePath = 'scaling_analysis.xlsx';
startRow = 1;
sheetComparisonsC28Star='LutealParametricComparisons';

human_conditions = conditions(contains(conditions, "human"));
%cs_avg_conditions = cell(1,length(conditions));
%c0_conditions = cell(1,length(conditions));
%m0_conditions = cell(1,length(conditions));
T_conditions_luteal = cell(1,length(conditions));
cs_avg_human_baselines_luteal = cell(1,length(human_conditions));

for i = 1:length(conditions)
    condition = conditions(i);
    [C, CI_avg_4C, CF_avg_4C, CE_avg_4C, CS_avg_4C, M_4C, N, Mass_remaining_4C, mR_inner, CB_4C, C_0, M_0, T, params] = solve_diffusion_5C(t, condition, drug, 'nonDimensional', true);
    
    startRow = writeExcelTable(T,filePath,sheetComparisonsC28Star,startRow,condition);
    %cs_avg_conditions{i} = CS_avg_4C;
    %c0_conditions{i} = C_0;
    %m0_conditions{i} = M_0;
    T_conditions_luteal{i} = T;

    % If the current condition is a human condition, store the cs_avg vector
    if contains(condition, "human")
        human_index = find(contains(human_conditions, condition));
        cs_avg_human_baselines_luteal{human_index} = CS_avg_4C;
    end

end

%% Parametric combinations - C28*/C0 calculation

% Extract human and macaque conditions directly
human_conditions = conditions(contains(conditions, "human"));
macaque_conditions = conditions(contains(conditions, "macaque"));

% Count the number of human and macaque conditions
numHumanConditions = length(human_conditions);
numMacaqueConditions = length(macaque_conditions);

% Identify human and macaque indices
human_indices = find(contains(conditions, "human"));
macaque_indices = find(contains(conditions, "macaque"));

% Calculate the total number of valid combinations
totalCombinations = numHumanConditions * numMacaqueConditions;

% Pre-allocate the comparisonResults cell array for efficiency
comparisonResults_luteal = cell(totalCombinations, 3); % 3 columns for 'Condition 1', 'Condition 2', and 'C0_optimal_Cstar28'
c0_optimal_cstar28_luteal = cell(numHumanConditions, numMacaqueConditions);
cs_avg_cstar28_luteal = cell(numHumanConditions, numMacaqueConditions);

% Initialize an index to fill in comparisonResults
resultIndex = 1;

% Iterate over each human and macaque condition
for i = 1:numHumanConditions
    for j = 1:numMacaqueConditions       
        % Get the indices for the current human and macaque conditions
        human_index = human_indices(i);
        macaque_index = macaque_indices(j);

        human_condition = conditions(human_index);
        macaque_condition = conditions(macaque_index);

        % Calculate the metric and run the model
        c0_optimal = T_conditions_luteal{human_index}.("C*28/C0")(4) / T_conditions_luteal{macaque_index}.("C*28/C0")(4);
        c0_optimal_cstar28_luteal{i,j} = c0_optimal;
        [~, ~, ~, ~, cstar28_optimal_cs_avg, ~, ~, ~, ~, ~, ~, ~, ~, ~] = solve_diffusion_5C(t, macaque_condition, drug, 'C_0', c0_optimal);
        cs_avg_cstar28_luteal{i,j} = cstar28_optimal_cs_avg;

        % Store the result in the pre-allocated cell array
        comparisonResults_luteal{resultIndex, 1} = human_condition;
        comparisonResults_luteal{resultIndex, 2} = macaque_condition;
        comparisonResults_luteal{resultIndex, 3} = c0_optimal_cstar28_luteal{i,j};

        % Increment the index
        resultIndex = resultIndex + 1;
    end
end

% Convert the results to a table
comparisonTable = cell2table(comparisonResults_luteal, 'VariableNames', {'Condition1', 'Condition2', 'c0_optimal_cstar28'});

% Write the table to the Excel sheet
startRow = writeExcelTable(comparisonTable, filePath, sheetComparisonsC28Star, startRow, 'C*28/C0 Comparisons');

%% RMS Minimization for Multiple Combinations

drug = "hydrophilic";
%conditions = ["human_avg", "human_large", "human_small", "macaque_avg", "macaque_large", "macaque_small"];
filePath = 'scaling_analysis.xlsx';
startRow = 1;
sheetMultipleParameterOptimization = 'LutealMultOptimizationMacaque';

% Identify human and macaque indices
human_indices = find(contains(conditions, "human"));
macaque_indices = find(contains(conditions, "macaque"));

% Get the number of human and macaque conditions
numHumanConditions = length(human_indices);
numMacaqueConditions = length(macaque_indices);

% Initialize cell arrays with the correct size
cs_avg_optimal_conditions_luteal = cell(numHumanConditions, numMacaqueConditions);
C0_rmse_logs_luteal = cell(numHumanConditions, numMacaqueConditions);
optimal_C0_values_luteal = cell(numHumanConditions, numMacaqueConditions);

% Iterate over all human and macaque combinations
for i = 1:numHumanConditions
    for j = 1:numMacaqueConditions
        % Get the indices for the current human and macaque conditions
        human_index = human_indices(i);
        macaque_index = macaque_indices(j);

        % Perform human calculation for reference
        condition_human = conditions(human_index);
        [~, ~, ~, ~, cs_avg_human, ~, ~, ~, ~, ~, ~, ~, ~] = solve_diffusion_5C(t, condition_human, drug, 'C_0', 1);

        % Set up optimization for macaque condition
        initialGuess_C0 = 1;
        error_func_args = {'time_pt', 28*24*3600};
        condition_macaque = conditions(macaque_index);

        % Run optimization for the current pair
        [optimal_C0, C0_rmse_log, cs_avg_macaque_optimal] = runOptimization(t, cs_avg_human, drug, condition_macaque, @rmse_custom, initialGuess_C0, error_func_args);

        % Store the results for plotting and analysis
        cs_avg_optimal_conditions_luteal{i, j} = cs_avg_macaque_optimal;
        C0_rmse_logs_luteal{i, j} = C0_rmse_log;
        optimal_C0_values_luteal{i, j} = [optimal_C0, C0_rmse_log(end, 2)]; % [C0_optimal, M0_optimal]

        % Create a table from the C0_rmse_log
        C0_rmse_table_luteal = array2table(C0_rmse_log, 'VariableNames', {'C0_macaque', 'M0_macaque', 'RMSE'});
        heading = sprintf("Optimal C0 for %s vs. %s: C0 = %.2f", condition_human, condition_macaque, optimal_C0);

        % Write the table to the Excel sheet
        startRow = writeExcelTable(C0_rmse_table_luteal, filePath, sheetMultipleParameterOptimization, startRow, heading);
    end
end

%% Writing Summary Results for RMSE optimization of multiple comparisons

% Preallocate the summaryResults cell array
summaryResults_luteal = cell(numHumanConditions * numMacaqueConditions, 4); % 4 columns for Condition1, Condition2, Optimal_C0, Optimal_M0

% Fill in the summaryResults with optimal values
resultIndex = 1;
for i = 1:numHumanConditions
    for j = 1:numMacaqueConditions
        % Store human and macaque conditions
        summaryResults_luteal{resultIndex, 1} = human_conditions{i};
        summaryResults_luteal{resultIndex, 2} = macaque_conditions{j};
        
        % Store the optimal C0 and M0 values
        summaryResults_luteal{resultIndex, 3} = optimal_C0_values_luteal{i, j}(1); % Optimal_C0
        summaryResults_luteal{resultIndex, 4} = optimal_C0_values_luteal{i, j}(2); % Optimal_M0
        
        % Increment the index
        resultIndex = resultIndex + 1;
    end
end

% Convert to table with appropriate column names
summaryTable_luteal = cell2table(summaryResults_luteal, 'VariableNames', {'Condition1', 'Condition2', 'Optimal_C0', 'Optimal_M0'});

% Write the summary table to the Excel sheet
startRow = writeExcelTable(summaryTable_luteal, filePath, sheetMultipleParameterOptimization, startRow, 'Summary of Optimal C0 and M0 for All Combinations - Luteal phase');

%% RUNNING ALLOMETRIC ESTIMATIONS

% Define the allometric C0 values matrix (numHumanConditions x numMacaqueConditions)
allometric_C0_values = [
    1.04, 1.87, 0.63;  % human_avg vs macaque_avg, macaque_large, macaque_small
    0.65, 1.17, 0.39;  % human_large vs macaque_avg, macaque_large, macaque_small
    1.96, 3.52, 1.17   % human_small vs macaque_avg, macaque_large, macaque_small
];

cs_avg_allometric_conditions = cell(numHumanConditions, numMacaqueConditions);

for i = 1:length(human_indices)
        human_index = human_indices(i);
    for j = 1:length(macaque_indices)
        macaque_index = macaque_indices(j);

        C0_allometric = allometric_C0_values(i, j);
        [~, ~, ~, ~, cs_avg_allometric, ~, ~, ~, ~, ~, ~, ~, ~, ~] = solve_diffusion_5C(t, conditions(macaque_index), drug, 'C_0', C0_allometric);
        cs_avg_allometric_conditions{i,j} = cs_avg_allometric;
    end
end


%% SAVING .MAT FILE

% Define the conditions array
%conditions = ["human_avg", "human_large", "human_small", "macaque_avg", "macaque_large", "macaque_small"];

% Create separate arrays for human and macaque conditions
human_conditions = conditions(contains(conditions, "human"));
macaque_conditions = conditions(contains(conditions, "macaque"));

% Save all relevant variables for plotting in a .mat file
outputFileName = 'pk_metrics_data_luteal.mat';
save(outputFileName, 't', 'human_conditions', 'macaque_conditions', ...
    'cs_avg_human_baselines_luteal', ... % HUMAN BASELINES 
    'c0_optimal_cstar28_luteal','cs_avg_cstar28_luteal', ... %C%28 ANALYSIS
    'optimal_C0_values_luteal', 'cs_avg_optimal_conditions_luteal', ... %RMSE ANALYSIS
    'allometric_C0_values', 'cs_avg_allometric_conditions', ... %ALLOMETRIC ANALYSIS
    'drug', '-v7.3');


%% PLOTTING RESULTS

% Assuming you have already created the variables:
% t: time vector
% conditions: array of condition names (e.g., ["human_avg", "human_large", "human_small", "macaque_avg", "macaque_large", "macaque_small"])
% cs_avg_conditions: cell array containing Cs/C0 results for each condition
% cs_avg_optimal_conditions: cell array containing RMSE minimization results for each human-macaque pair
% allometric_C0_values: matrix of allometric C0 values (3x3 matrix for 3 human and 3 macaque conditions)
% solve_diffusion_5C: function handle to your model function
% drug: string representing the drug type (e.g., "hydrophilic")

plotPKComparisons(t, conditions, cs_avg_human_baselines_luteal, cs_avg_cstar28_luteal, cs_avg_optimal_conditions_luteal, cs_avg_allometric_conditions);
%plotPKVariability(t, conditions, cs_avg_human_baselines_luteal, cs_avg_cstar28_luteal, cs_avg_optimal_conditions_luteal, cs_avg_allometric_conditions)

%% PLOTTING PK vs R

drug = "hydrophilic";
conditions =["human_avg_luteal","human_large_luteal","macaque_avg_luteal","macaque_small_luteal"];
filePath = 'scaling_analysis.xlsx';

human_conditions = conditions(contains(conditions, "human"));
%cs_avg_conditions = cell(1,length(conditions));
%c0_conditions = cell(1,length(conditions));
%m0_conditions = cell(1,length(conditions));
T_conditions_luteal = cell(1,length(conditions));
cs_avg_human_baselines_luteal = cell(1,length(human_conditions));

for i = 1:length(conditions)
    condition = conditions(i);

    % If the current condition is a human condition, store the cs_avg vector
    if contains(condition, "human")
        [C, CI_avg_4C, CF_avg_4C, CE_avg_4C, CS_avg_4C, M_4C, N, Mass_remaining_4C, mR_inner, CB_4C, C_0, M_0, T, params] = solve_diffusion_5C(t, condition, drug, 'nonDimensional', true);
    end
    max(C(end,:))

end

%% STEP ANALYSIS

% Parametric combinations - Running all conditions non-dimensionalized

drug = "hydrophilic";
isNonDimensionalized = true;
conditions =["human_step_increasing_size","macaque_step_increasing_size"];
filePath = 'scaling_analysis.xlsx';
startRow = 1;
sheetComparisonsC28Star='step_analysis_c28star';
numSteps = 20;

human_conditions = conditions(contains(conditions, "human"));
T_conditions = cell(length(conditions), numSteps);  % Store T for each step
cs_avg_human_baselines = cell(length(human_conditions), numSteps);  % Store baseline for each step


for j = 1:numSteps
    for i = 1:length(conditions)
        condition = conditions(i);
        [C, CI_avg_4C, CF_avg_4C, CE_avg_4C, CS_avg_4C, M_4C, N, Mass_remaining_4C, mR_inner, CB_4C, C_0, M_0, T, params] = solve_diffusion_5C(t, condition, drug, 'nonDimensional', true, 'numSteps', numSteps, 'currentStep', j);
        
        startRow = writeExcelTable(T,filePath,sheetComparisonsC28Star,startRow,sprintf("%s - step_%d", condition, j));
        startRow = writeExcelTable(params,filePath,sheetComparisonsC28Star,startRow,sprintf("%s - step_%d, parameters",condition, j));
        T_conditions{i,j} = T;
    
        % If the current condition is a human condition, store the cs_avg vector
        if contains(condition, "human")
            human_index = find(contains(human_conditions, condition));
            cs_avg_human_baselines{human_index, j} = CS_avg_4C;
        end
    end
end

%% Parametric Combinations - C28*/C0 Calculation


% Extract human and macaque conditions directly
human_conditions = conditions(contains(conditions, "human"));
macaque_conditions = conditions(contains(conditions, "macaque"));

% Count the number of human and macaque conditions
numHumanConditions = length(human_conditions);
numMacaqueConditions = length(macaque_conditions);

% Identify human and macaque indices
human_indices = find(contains(conditions, "human"));
macaque_indices = find(contains(conditions, "macaque"));

% Pre-allocate the comparisonResults cell array to store results
% Columns: Human Condition, Macaque Condition, Step for Human, Step for Macaque, Optimal C0 Ratio
comparisonResults = cell(numHumanConditions * numMacaqueConditions * numSteps * numSteps, 5);
cs_avg_cstar28 = cell(numHumanConditions, numMacaqueConditions, numSteps, numSteps);           % Store cs_avg for each combination
c0_optimal_cstar28 = cell(numHumanConditions, numMacaqueConditions, numSteps, numSteps);       % Store c0_optimal for each combination

% Initialize result index
resultIndex = 1;

% Iterate through each human condition
for i = 1:numHumanConditions
    
    % Iterate through each macaque condition
    for j = 1:numMacaqueConditions
        % human_condition = string(human_conditions{humanIndex});
        % macaque_condition = string(macaque_conditions{macaqueIndex});
        human_index = human_indices(i);
        macaque_index = macaque_indices(j);
        human_condition = conditions(human_index);
        macaque_condition = conditions(macaque_index);

        % Nested loops to cover all combinations of steps for human and macaque conditions
        for human_step = 1:numSteps
            for macaque_step = 1:numSteps
                % Retrieve the specific "C*28/C0" values at each step
                human_Cstar28 = T_conditions{human_index, human_step}.("C*28/C0")(4);
                macaque_Cstar28 = T_conditions{macaque_index, macaque_step}.("C*28/C0")(4);

                % Calculate the optimal C0 ratio for the current combination
                c0_optimal = human_Cstar28 / macaque_Cstar28;  % Get optimal C0 ratio

                % Store the optimal C0 value in the array
                c0_optimal_cstar28{human_index, macaque_index, human_step, macaque_step} = c0_optimal;

                % Run the model with this optimal C0 value for the macaque condition at the given step
                [~, ~, ~, ~, cs_avg_optimal, ~, ~, ~, ~, ~, ~, ~, ~] = solve_diffusion_5C( ...
                    t, macaque_condition, drug, 'C_0', c0_optimal, 'numSteps', numSteps, 'currentStep', macaque_step);

                % Store the resulting cs_avg in the array
                cs_avg_cstar28{human_index, macaque_index, human_step, macaque_step} = cs_avg_optimal;

                % Store results in comparisonResults
                comparisonResults{resultIndex, 1} = human_condition;        % Human condition name
                comparisonResults{resultIndex, 2} = macaque_condition;      % Macaque condition name
                comparisonResults{resultIndex, 3} = human_step;              % Step number for human condition
                comparisonResults{resultIndex, 4} = macaque_step;            % Step number for macaque condition
                comparisonResults{resultIndex, 5} = c0_optimal;             % Optimal C0 ratio (C*28/C0)

                % Increment the result index for the next entry
                resultIndex = resultIndex + 1;
            end
        end
    end
end

% Convert the results into a table with descriptive column names
comparisonTable = cell2table(comparisonResults, ...
    'VariableNames', {'HumanCondition', 'MacaqueCondition', 'Step_Human', 'Step_Macaque', 'C0_optimal_Cstar28'});

%% Write .mat file
% Write the table to the Excel sheet
outputFileName = 'pk_metrics_data_step_increasing_size.mat';
save(outputFileName, 't', 'human_conditions', 'macaque_conditions', ...
    'cs_avg_human_baselines', ... % HUMAN BASELINES 
    'c0_optimal_cstar28','cs_avg_cstar28', ... % C*28 ANALYSIS
    'comparisonResults', ...
    'drug', '-v7.3');

%% plot histogram

color = [255/255,127/255,80/255];
plotHistogramOptimalC0(c0_optimal_cstar28,'fitType',"Lognormal",'ylabel',"C*28/C0 PDF",'color',color);
%plotHistogramOptimalC0(c0_optimal_cstar28,'fitType',"Gamma");
%plotHistogramOptimalC0(c0_optimal_cstar28,'fitType',"Weibull");

%% allometric equivalent analysis

% Define constants and conditions
Vrh = 6.61;  % Human reference volume
Vrm = 1.69;  % Macaque reference volume
numSteps = 20;  % Number of step-increasing sizes

% Define the step-increasing sizes for human and macaque conditions
A_F_human = linspace(80, 120, numSteps);
h_S_human = linspace(0.1, 0.2, numSteps);
A_F_macaque = linspace(32, 48, numSteps);
h_S_macaque = linspace(0.075, 0.15, numSteps);

% Calculate tissue volumes
Vsh = A_F_human .* h_S_human;  % Human tissue volumes for each step
Vsm = A_F_macaque .* h_S_macaque;  % Macaque tissue volumes for each step

% Compute the allometric C0 values matrix
allometric_C0_values = (Vrh ./ Vsh') ./ (Vrm ./ Vsm);  % Result is numSteps x numSteps

%% Write .mat file
% Write the table to the Excel sheet
outputFileName = 'pk_metrics_data_step_increasing_size.mat';
save(outputFileName, 't', 'human_conditions', 'macaque_conditions', ...
    'cs_avg_human_baselines', ... % HUMAN BASELINES 
    'c0_optimal_cstar28','cs_avg_cstar28', ... % C*28 ANALYSIS
    'allometric_C0_values', ... % ALLOMETRIC ANALYSIS
    'comparisonResults', ...
    'drug', '-v7.3');

%% plot allometric histogram
color = [80/255,127/255,255/255];
plotHistogramOptimalC0(num2cell(allometric_C0_values),'fitType',"Lognormal",'ylabel',"Allometric PDF", 'color', color);

%% 60mg M0 analysis

% Parametric combinations - Running all conditions non-dimensionalized

t = linspace(0, 60*60*24*144, 10000);
t_day = t./(60*60*24);
drug = "hydrophilic";
isNonDimensionalized = true;
conditions =["human_avg","human_large","human_small","macaque_avg","macaque_large","macaque_small"];
filePath = 'scaling_analysis.xlsx';
startRow = 1;
sheetComparisonsC28Star_60mg='HumanVsMacaque60mgLoading';

T_conditions = cell(1,length(conditions));
cs_avg_all_conditions_60load = cell(1,length(conditions));

for i = 1:length(conditions)
    condition = conditions(i);
    [C, CI_avg_4C, CF_avg_4C, CE_avg_4C, CS_avg_4C, M_4C, N, Mass_remaining_4C, mR_inner, CB_4C, C_0, M_0, T, params] = solve_diffusion_5C(t, condition, drug, 'M_0', 60);
    
    startRow = writeExcelTable(T,filePath,sheetComparisonsC28Star_60mg,startRow,condition);
    T_conditions{i} = T;
    cs_avg_all_conditions_60load{i} = CS_avg_4C/C_0;

end

%% PLOTTING COMPARISONS

plotUnoptimizedPKComparisons(t,conditions,cs_avg_all_conditions_60load);

%% STEP ANALYSIS - CYCLE PHASE SHEEP, MACAQUE, HUMAN

% Parametric combinations - Running all conditions non-dimensionalized
t = linspace(0, 60*60*24*144, 10000);
t_day = t./(60*60*24);
drug = "hydrophilic";
isNonDimensionalized = true;
conditions = ["human_step_increasing_size_follicular", "human_step_increasing_size_midcycle", "human_step_increasing_size_luteal", ...
              "macaque_step_increasing_size_follicular", "macaque_step_increasing_size_midcycle", "macaque_step_increasing_size_luteal", ...
              "sheep_step_increasing_size_follicular", "sheep_step_increasing_size_midcycle", "sheep_step_increasing_size_luteal"];
filePath = 'scaling_analysis.xlsx';
startRow = 1;
sheetComparisonsC28Star = 'step_c28star_cycle_size';
numSteps = 2; % Adjust as needed for debugging or full analysis

% Extract conditions for each organism
human_conditions = conditions(contains(conditions, "human"));
macaque_conditions = conditions(contains(conditions, "macaque"));
sheep_conditions = conditions(contains(conditions, "sheep"));

% Pre-allocate cells to store results
T_conditions = cell(length(conditions), numSteps); % Store T for each condition and step
cs_avg_results = cell(length(conditions), numSteps); % Store baseline for each condition and step

% Run the diffusion model for each condition and step
for step = 1:numSteps
    for i = 1:length(conditions)
        condition = conditions(i);
        [C, CI_avg_4C, CF_avg_4C, CE_avg_4C, CS_avg_4C, M_4C, N, Mass_remaining_4C, mR_inner, CB_4C, C_0, M_0, T, params] = solve_diffusion_5C(t, condition, drug, 'nonDimensional', true, 'numSteps', numSteps, 'currentStep', step);
        startRow = writeExcelTable(T, filePath, sheetComparisonsC28Star, startRow, sprintf("%s - step_%d", condition, step));
        startRow = writeExcelTable(params, filePath, sheetComparisonsC28Star, startRow, sprintf("%s - step_%d, parameters", condition, step));
        T_conditions{i, step} = T;
        cs_avg_results{i, step} = CS_avg_4C;
    end
end

%% PARAMETRIC COMBINATIONS - C28*/C0 Calculation

% Precompute indices for each condition
human_indices = find(contains(conditions, "human"));
macaque_indices = find(contains(conditions, "macaque"));
sheep_indices = find(contains(conditions, "sheep"));

% Create a grid for step combinations
[human_steps, macaque_steps, sheep_steps] = ndgrid(1:numSteps, 1:numSteps, 1:numSteps);

% Flatten step grids into vectors for vectorized processing
human_steps = human_steps(:);
macaque_steps = macaque_steps(:);
sheep_steps = sheep_steps(:);

% Pre-allocate results
numCombinations = numel(human_steps) * numHumanConditions * numMacaqueConditions * numSheepConditions;
comparisonResults = cell(numCombinations, 8);
cs_avg_cstar28 = cell(numCombinations, 1);
c0_optimal_cstar28 = cell(numCombinations, 1);

% Initialize result index
resultIndex = 1;

% Iterate through all combinations of human, macaque, and sheep conditions
for humanIdx = 1:numHumanConditions
    for macaqueIdx = 1:numMacaqueConditions
        for sheepIdx = 1:numSheepConditions
            human_condition = conditions(human_indices(humanIdx));
            macaque_condition = conditions(macaque_indices(macaqueIdx));
            sheep_condition = conditions(sheep_indices(sheepIdx));

            % Retrieve the "C*28/C0" values for all steps in one go
            human_Cstar28 = arrayfun(@(step) T_conditions{human_indices(humanIdx), step}.("C*28/C0")(4), human_steps);
            macaque_Cstar28 = arrayfun(@(step) T_conditions{macaque_indices(macaqueIdx), step}.("C*28/C0")(4), macaque_steps);
            sheep_Cstar28 = arrayfun(@(step) T_conditions{sheep_indices(sheepIdx), step}.("C*28/C0")(4), sheep_steps);

            % Calculate the optimal C0 ratios
            c0_human_macaque = human_Cstar28 ./ macaque_Cstar28;
            c0_human_sheep = human_Cstar28 ./ sheep_Cstar28;

            % Run the model for each combination
            for stepIdx = 1:numel(human_steps)
                cs_avg_optimal = solve_diffusion_5C( ...
                    t, macaque_condition, drug, 'C_0', c0_human_macaque(stepIdx), ...
                    'numSteps', numSteps, 'currentStep', macaque_steps(stepIdx));
                cs_avg_cstar28{resultIndex} = cs_avg_optimal;

                % Store results in comparisonResults
                comparisonResults{resultIndex, 1} = human_condition;
                comparisonResults{resultIndex, 2} = macaque_condition;
                comparisonResults{resultIndex, 3} = sheep_condition;
                comparisonResults{resultIndex, 4} = human_steps(stepIdx);
                comparisonResults{resultIndex, 5} = macaque_steps(stepIdx);
                comparisonResults{resultIndex, 6} = sheep_steps(stepIdx);
                comparisonResults{resultIndex, 7} = c0_human_macaque(stepIdx);
                comparisonResults{resultIndex, 8} = c0_human_sheep(stepIdx);

                % Increment result index
                resultIndex = resultIndex + 1;
            end
        end
    end
end

% Convert results into a table with descriptive column names
comparisonTable = cell2table(comparisonResults, ...
    'VariableNames', {'HumanCondition', 'MacaqueCondition', 'SheepCondition', 'Step_Human', 'Step_Macaque', 'Step_Sheep', ...
                      'C0_optimal_Human_Macaque', 'C0_optimal_Human_Sheep'});

% Write the table to the Excel sheet
startRow = writeExcelTable(comparisonTable, filePath, sheetComparisonsC28Star, startRow, 'Step Analysis - All Conditions');


%% Write .mat file
% Write the table to the Excel sheet
outputFileName = 'pk_metrics_data_step_increasing_size_phase.mat';
save(outputFileName, 't', 'human_conditions', 'macaque_conditions', 'sheep_conditions', ...
    'cs_avg_human_baselines', ... % HUMAN BASELINES 
    'c0_optimal_cstar28','cs_avg_cstar28', ... % C*28 ANALYSIS
    'comparisonResults', ...
    'drug', '-v7.3');

%% plot histogram

color = [255/255,127/255,80/255];
plotHistogramOptimalC0(c0_optimal_cstar28,'fitType',"Lognormal",'ylabel',"C*28/C0 PDF",'color',color);
%plotHistogramOptimalC0(c0_optimal_cstar28,'fitType',"Gamma");
%plotHistogramOptimalC0(c0_optimal_cstar28,'fitType',"Weibull");

%% allometric equivalent analysis

% Define constants and conditions
Vrh = 6.61;  % Human reference volume
Vrm = 1.69;  % Macaque reference volume
numSteps = 20;  % Number of step-increasing sizes

% Define the step-increasing sizes for human and macaque conditions
A_F_human = linspace(80, 120, numSteps);
h_S_human = linspace(0.1, 0.2, numSteps);
A_F_macaque = linspace(32, 48, numSteps);
h_S_macaque = linspace(0.075, 0.15, numSteps);

% Calculate tissue volumes
Vsh = A_F_human .* h_S_human;  % Human tissue volumes for each step
Vsm = A_F_macaque .* h_S_macaque;  % Macaque tissue volumes for each step

% Compute the allometric C0 values matrix
allometric_C0_values = (Vrh ./ Vsh') ./ (Vrm ./ Vsm);  % Result is numSteps x numSteps

%% Write .mat file
% Write the table to the Excel sheet
outputFileName = 'pk_metrics_data_step_increasing_size.mat';
save(outputFileName, 't', 'human_conditions', 'macaque_conditions', ...
    'cs_avg_human_baselines', ... % HUMAN BASELINES 
    'c0_optimal_cstar28','cs_avg_cstar28', ... % C*28 ANALYSIS
    'allometric_C0_values', ... % ALLOMETRIC ANALYSIS
    'comparisonResults', ...
    'drug', '-v7.3');

%% plot allometric histogram
color = [80/255,127/255,255/255];
plotHistogramOptimalC0(num2cell(allometric_C0_values),'fitType',"Lognormal",'ylabel',"Allometric PDF", 'color', color);

