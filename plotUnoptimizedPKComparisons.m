function plotUnoptimizedPKComparisons(t, conditions, cs_avgs)
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
    colors = lines(2); % Generate distinct colors for 4 lines

    % Plot human_avg line and variability shading
    plot(t_days, cs_avgs{human_avg_index}, 'Color', colors(1, :), 'LineWidth', 4, 'DisplayName', 'Human Avg');
    fill_between(t_days, cs_avgs{human_small_index}, cs_avgs{human_large_index}, colors(1, :), 0.1);

    % Plot macaque_avg line and variability shading
    plot(t_days, cs_avgs{macaque_avg_index}, 'Color', colors(2, :), 'LineWidth', 4, 'DisplayName', 'Macaque Avg');
    fill_between(t_days, cs_avgs{macaque_small_index}, cs_avgs{macaque_large_index}, colors(2, :), 0.3);

    legendLabels={strcat("Human"),'',strcat("Macaque"),''};
    % Customize the plot
    set(gca, 'FontSize', 28);
    xlabel('Time (days)', 'FontSize', 36);
    ylabel('Cs/C0', 'FontSize', 36);
    xlim([0, 28]);
    legend(legendLabels,'Location', 'northeast');
    title('60mg IVR Loading');
    hold off;
end

function fill_between(x, y1, y2, color, alpha)
    % Fills the area between two curves (y1 and y2) with a transparent color
    fill([x, fliplr(x)], [y1, fliplr(y2)], color, 'FaceAlpha', alpha, 'EdgeColor', 'none');
end
