function [optimal_C0, C0_rmse_log, cs_avg_experiment_optimal] = runOptimization(t, cs_avg_reference, drug, condition, error_func, initialGuess, error_func_args)
    % Initialize a local variable to log C_0 and RMSE values
    C0_rmse_log = [];
    cs_avg_experiment_optimal = []; % Variable to store the optimal cs_avg_experiment
    best_rmse = Inf; % Initial value to track the best RMSE


    % Create a function handle that logs data
    objectiveFunction = @(C_0_experiment) logObjective(C_0_experiment);

    % Nested function to log and calculate RMSE
    function rmseValue = logObjective(C_0_experiment)
        if C_0_experiment <= 0
            rmseValue = 1e6; % Assign a large penalty value for invalid C_0
            return;
        end

        [~, ~, ~, ~, cs_avg_experiment, ~, ~, ~, ~, ~, ~, M_0_experiment, ~] = solve_diffusion_5C(t, condition, drug, 'C_0', C_0_experiment);
        

        rmseValue = error_func(t, cs_avg_reference, cs_avg_experiment, error_func_args{:});
        
        % Append to the log
        C0_rmse_log = [C0_rmse_log; C_0_experiment, M_0_experiment, rmseValue];

        % Update cs_avg_experiment_optimal if the current RMSE is better than the best so far
        if rmseValue < best_rmse
            best_rmse = rmseValue;
            cs_avg_experiment_optimal = cs_avg_experiment;
        end
    end

    % Run fminsearch to find the optimal C_0 for the macaque
    optimal_C0 = fminsearch(objectiveFunction, initialGuess);
end