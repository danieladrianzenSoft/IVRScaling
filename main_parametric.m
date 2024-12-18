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

% Define the number of conditions for each organism
numHumanConditions = length(human_conditions);
numMacaqueConditions = length(macaque_conditions);
numSheepConditions = length(sheep_conditions);

% Pre-allocate cells to store results
T_conditions = cell(length(conditions), numSteps); % Store T for each condition and step

% Run the diffusion model for each condition and step
for step = 1:numSteps
    for i = 1:length(conditions)
        condition = conditions(i);
        [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, T, params] = solve_diffusion_5C(t, condition, drug, ...
            'nonDimensional', true, 'numSteps', numSteps, 'currentStep', step);
        startRow = writeExcelTable(T, filePath, sheetComparisonsC28Star, startRow, sprintf("%s - step_%d", condition, step));
        startRow = writeExcelTable(params, filePath, sheetComparisonsC28Star, startRow, sprintf("%s - step_%d, parameters", condition, step));
        T_conditions{i, step} = T;
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
comparisonResults = cell(numCombinations, 9);

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
            c0_macaque_sheep = macaque_Cstar28 ./ sheep_Cstar28;

            % Store results in comparisonResults
            for stepIdx = 1:numel(human_steps)
                comparisonResults{resultIndex, 1} = human_condition;
                comparisonResults{resultIndex, 2} = macaque_condition;
                comparisonResults{resultIndex, 3} = sheep_condition;
                comparisonResults{resultIndex, 4} = human_steps(stepIdx);
                comparisonResults{resultIndex, 5} = macaque_steps(stepIdx);
                comparisonResults{resultIndex, 6} = sheep_steps(stepIdx);
                comparisonResults{resultIndex, 7} = c0_human_macaque(stepIdx);
                comparisonResults{resultIndex, 8} = c0_human_sheep(stepIdx);
                comparisonResults{resultIndex, 9} = c0_macaque_sheep(stepIdx); % Add new ratio
                resultIndex = resultIndex + 1;
            end
        end
    end
end

% Convert results into a table with descriptive column names
comparisonTable = cell2table(comparisonResults, ...
    'VariableNames', {'HumanCondition', 'MacaqueCondition', 'SheepCondition', 'Step_Human', 'Step_Macaque', 'Step_Sheep', ...
                      'C0_optimal_Human_Macaque', 'C0_optimal_Human_Sheep', 'C0_optimal_Macaque_Sheep'});

% Write the table to the Excel sheet
startRow = writeExcelTable(comparisonTable, filePath, sheetComparisonsC28Star, startRow, 'Step Analysis - All Conditions');

%% ALLOMETRIC EQUIVALENT ANALYSIS

% Define constants and conditions
Vrh = 6.61;  % Human reference volume
Vrm = 1.69;  % Macaque reference volume
Vrs = 6.61;  % Sheep reference volume

% Define the step-increasing sizes for human and macaque conditions
A_F_human = linspace(70, 130, numSteps);
h_S_human = linspace(.126, .154, numSteps);
A_F_macaque = linspace(32, 59, numSteps);
h_S_macaque = linspace(0.09, 0.11, numSteps);
A_F_sheep = linspace(53, 97, numSteps);
h_S_sheep = linspace(.117, .143, numSteps);

% Calculate tissue volumes
Vsh = A_F_human .* h_S_human;  % Human tissue volumes for each step
Vsm = A_F_macaque .* h_S_macaque;  % Macaque tissue volumes for each step
Vss = A_F_sheep .* h_S_sheep;  % Macaque tissue volumes for each step

% Compute the allometric C0 values matrix
allometric_C0_Human_Macaque = (Vrh ./ Vsh') ./ (Vrm ./ Vsm);  % Result is numSteps x numSteps
allometric_C0_Human_Sheep = (Vrh ./ Vsh') ./ (Vrs ./ Vss);  % Result is numSteps x numSteps
allometric_C0_Macaque_Sheep = (Vrm ./ Vsm') ./ (Vrs ./ Vss);  % Result is numSteps x numSteps

%% Write .mat file
% Write the table to the Excel sheet
outputFileName = 'pk_metrics_data_step_increasing_size_phase.mat';
save(outputFileName, 't', 'human_conditions', 'macaque_conditions', 'sheep_conditions', ... % CONDITIONS
    'allometric_C0_Human_Macaque', 'allometric_C0_Human_Sheep', 'allometric_C0_Macaque_Sheep', ... % ALLOMETRIC ANALYSIS
    'comparisonResults', ... % COMPARISON RESULTS
    'drug', '-v7.3');


