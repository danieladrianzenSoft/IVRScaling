function plotPKComparisons(t, conditions, cs_avg_human_baselines, cs_avg_cstar28, cs_avg_optimal_conditions, cs_avg_allometric_conditions)
    % Identify human and macaque indices
    human_indices = find(contains(conditions, "human"));
    macaque_indices = find(contains(conditions, "macaque"));
    %cs_avg_optimal_conditions_sub = cs_avg_optimal_conditions(1:3, 4:6);
    cs_avg_optimal_conditions_sub = cs_avg_optimal_conditions;
    
    t_days = t / (24 * 3600); % Convert time to days

    % Iterate over all human and macaque combinations for plotting
    for i = 1:length(human_indices)
        human_index = human_indices(i);
        
        for j = 1:length(macaque_indices)
            macaque_index = macaque_indices(j);

            % Create a new figure for each combination
            figure();
            hold on;
            
            % Plot the human baseline
            plot(t_days, cs_avg_human_baselines{human_index}, 'k--', 'LineWidth', 4);
            legendLabels = {strcat("Human Baseline")};
            
            % Plot C*28/C0 results
            plot(t_days, cs_avg_cstar28{i, j}, 'LineWidth', 4);
            legendLabels{end+1} = strcat("Rule 1 - C*28/C0");

            % Plot RMSE minimization results if available
            if ~isempty(cs_avg_optimal_conditions_sub{i, j})
                plot(t_days, cs_avg_optimal_conditions_sub{i, j}, 'LineWidth', 4);
                legendLabels{end+1} = strcat("Rule 2 - RMSE minimization");
            end
            
            % Plot allometric scaling results using the provided C0 values
            %C0_allometric = allometric_C0_values(i, j);
            %[~, ~, ~, ~, cs_avg_allometric, ~, ~, ~, ~, ~, ~, ~, ~] = solve_diffusion_5C(t, conditions(macaque_index), drug, 'C_0', C0_allometric);
            plot(t_days, cs_avg_allometric_conditions{i, j}, 'LineWidth', 4, 'DisplayName', ['Rule 3 - {' conditions{macaque_index} ' - Allometric}']);
            legendLabels{end+1} = strcat('Rule 3 - Allometric scaling');

            % Customize the plot
            set(gca, 'FontSize', 28);
            xlabel('Time (days)', 'FontSize', 36);
            ylabel('Cs/C0', 'FontSize', 36);
            xlim([0,28]);
            legend(legendLabels, 'Interpreter', 'none', 'Location', 'northeast');
            title(sprintf('%s vs. %s', conditions(human_index), conditions(macaque_index)), 'Interpreter', 'none');
            hold off;
        end
    end
end
