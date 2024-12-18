function plotPKVariability(t, conditions, cs_avg_human_baselines, cs_avg_cstar28, cs_avg_optimal_conditions, cs_avg_allometric_conditions)
    % Identify human and macaque indices
    human_indices = find(contains(conditions, "human"));
    macaque_indices = find(contains(conditions, "macaque"));
    
    % Get the number of human and macaque conditions
    numHumanConditions = length(human_indices);
    numMacaqueConditions = length(macaque_indices);

    % Identify indices for human and macaque conditions
    human_avg_index = find(contains(conditions, "human_avg"));
    human_small_index = find(contains(conditions, "human_small"));
    human_large_index = find(contains(conditions, "human_large"));
    
    macaque_avg_index = find(contains(conditions, "macaque_avg"));
    macaque_small_index = find(contains(conditions, "macaque_small"));
    macaque_large_index = find(contains(conditions, "macaque_large"));
    
    % Ensure correct indexing in the results matrices
    t_days = t / (24 * 3600); % Convert time to days

    % Create the figure
    figure();
    hold on;

    % Define colors for the lines and shaded regions
    colors = lines(4); % Generate distinct colors for 4 lines

    % Plot human_avg line and variability shading
    plot(t_days, cs_avg_human_baselines{human_avg_index}, 'Color', colors(1, :), 'LineWidth', 2, 'DisplayName', 'Human Avg');
    fill_between(t_days, cs_avg_human_baselines{human_small_index}, cs_avg_human_baselines{human_large_index}, colors(1, :), 0.1);

    % Plot macaque_avg C*28/C0 line and variability shading
    plot(t_days, cs_avg_cstar28{1, 1}, 'Color', colors(2, :), 'LineWidth', 2, 'DisplayName', 'Macaque Avg - C*28/C0');
    fill_between(t_days, cs_avg_cstar28{human_avg_index, macaque_small_index - numHumanConditions}, cs_avg_cstar28{human_avg_index, macaque_large_index - numHumanConditions}, colors(2, :), 0.3);

    % Plot macaque_avg RMSE line and variability shading
    plot(t_days, cs_avg_optimal_conditions{1, 1}, 'Color', colors(3, :), 'LineWidth', 2, 'DisplayName', 'Macaque Avg - RMSE');
    fill_between(t_days, cs_avg_optimal_conditions{human_avg_index, macaque_small_index - numHumanConditions}, cs_avg_optimal_conditions{human_avg_index, macaque_large_index - numHumanConditions}, colors(3, :), 0.3);

    % Plot macaque_avg Allometric line and variability shading
    plot(t_days, cs_avg_allometric_conditions{1, 1}, 'Color', colors(4, :), 'LineWidth', 2, 'DisplayName', 'Macaque Avg - Allometric');
    fill_between(t_days, cs_avg_allometric_conditions{human_avg_index, macaque_small_index - numHumanConditions}, cs_avg_allometric_conditions{human_avg_index, macaque_large_index - numHumanConditions}, colors(4, :), 0.3);

    legendLabels={strcat("Human"),'',strcat("Macaque - C*28/C0 optimization"),'',strcat("Macaque - RMSE minimization"),'',strcat('Macaque - Allometric scaling')};
    % Customize the plot
    set(gca, 'FontSize', 28);
    xlabel('Time (days)', 'FontSize', 36);
    ylabel('Cs/C0', 'FontSize', 36);
    xlim([0, 28]);
    legend(legendLabels,'Location', 'northeast');
    title('PK Metrics with Variability Bands');
    hold off;
end

function fill_between(x, y1, y2, color, alpha)
    % Fills the area between two curves (y1 and y2) with a transparent color
    fill([x, fliplr(x)], [y1, fliplr(y2)], color, 'FaceAlpha', alpha, 'EdgeColor', 'none');
end
