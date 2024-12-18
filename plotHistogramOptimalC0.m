function plotHistogramOptimalC0(c0_optimal_values, varargin)
    p = inputParser;
    p.addRequired("c0_optimal_values", @(x) length(x)>=1);
    p.addOptional("fitType", "Lognormal", @isstring);
    p.addOptional("title", "", @isstring);
    p.addOptional("ylabel", "PDF", @isstring);
    p.addOptional("color", [255/255,127/255,80/255], @(x) length(x)==3);
    p.parse(c0_optimal_values,varargin{:});
    fitType = p.Results.fitType;
    plotTitle = p.Results.title;
    ylabelText = p.Results.ylabel;
    color = p.Results.color;

    % Flatten the cell array to extract all non-empty optimal C0 values
    all_c0_optimal = cell2mat(c0_optimal_values(:));  % Convert cell array to matrix
    
    % Remove any empty or NaN entries, if they exist
    all_c0_optimal = all_c0_optimal(~isnan(all_c0_optimal) & all_c0_optimal ~= 0);

    % Create the histogram
    figure;
    h = histogram(all_c0_optimal, 'Normalization', 'pdf', 'FaceColor', color, 'FaceAlpha', 0.3);
    h.EdgeColor = 'none';
    xlabel('Optimal C_{0,m} / C_{0,h}', 'FontSize', 36);
    ylabel(ylabelText, 'FontSize', 36);
    set(gca, 'FontSize', 28);
    grid on;
    ylim([0,1.5])
    xlim([0,3.5])

    % Optional: Show basic statistics
    meanC0 = mean(all_c0_optimal)
    medianC0 = median(all_c0_optimal)
    stdC0 = std(all_c0_optimal);
    hold on;
    xline(meanC0, 'k', 'LineWidth', 4, 'DisplayName', sprintf('Mean = %.2f', meanC0));
    % xline(medianC0, 'color', [0.3, 0.3, 0.8], 'LineWidth', 2, 'DisplayName', sprintf('Median = %.2f', medianC0));

    % Generate x-values for PDF calculation
    x_values = linspace(min(all_c0_optimal), max(all_c0_optimal), 100);

    % Fit the distributions
    lognormal_pd = fitdist(all_c0_optimal, 'Lognormal');
    gamma_pd = fitdist(all_c0_optimal, 'Gamma');
    weibull_pd = fitdist(all_c0_optimal, 'Weibull');

    % Calculate and plot 95% CI for the Lognormal distribution mean
    mean_lognormal = exp(lognormal_pd.mu + (lognormal_pd.sigma^2) / 2);
    std_lognormal = sqrt((exp(lognormal_pd.sigma^2) - 1) * exp(2 * lognormal_pd.mu + lognormal_pd.sigma^2));
    se_lognormal = std_lognormal / sqrt(length(all_c0_optimal));
    ci_lognormal = [mean_lognormal - 1.96 * se_lognormal, mean_lognormal + 1.96 * se_lognormal];

    % Calculate and plot 95% CI for the Gamma distribution mean
    mean_gamma = gamma_pd.a * gamma_pd.b;
    std_gamma = sqrt(gamma_pd.a * (gamma_pd.b ^ 2));
    se_gamma = std_gamma / sqrt(length(all_c0_optimal));
    ci_gamma = [mean_gamma - 1.96 * se_gamma, mean_gamma + 1.96 * se_gamma];

    % Calculate and plot 95% CI for the Weibull distribution mean
    mean_weibull = weibull_pd.B * gamma(1 + 1 / weibull_pd.A);
    var_weibull = (weibull_pd.B ^ 2) * (gamma(1 + 2 / weibull_pd.A) - (gamma(1 + 1 / weibull_pd.A))^2);
    std_weibull = sqrt(var_weibull);
    se_weibull = std_weibull / sqrt(length(all_c0_optimal));
    ci_weibull = [mean_weibull - 1.96 * se_weibull, mean_weibull + 1.96 * se_weibull];
    
    y_limits = ylim; % Gets the current y-axis limits of the plot
    if strcmp(fitType, "Gamma") == 1
        plot(x_values, pdf(gamma_pd, x_values), 'color', color, 'LineWidth', 4);
        pd = gamma_pd;
        legendLabels = {'','Mean', 'Gamma fit', '2nd Quartile', '3rd Quartile'};
        %plot([ci_gamma(1), ci_gamma(1)], y_limits, '--k', 'LineWidth', 2, 'DisplayName', 'Lognormal 95% CI');
        %plot([ci_gamma(2), ci_gamma(2)], y_limits, '--k', 'LineWidth', 2);
        %text(ci_gamma(1), max(ylim) * 0.1, sprintf('%.2f', ci_gamma(1)), 'FontSize', 20, ...
        %    'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', 'Color', 'k');
        %text(ci_gamma(2), max(ylim) * 0.1, sprintf('%.2f', ci_gamma(2)), 'FontSize', 20, ...
        %    'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'Color', 'k');
    elseif strcmp(fitType, "Weibull") == 1
        plot(x_values, pdf(weibull_pd, x_values), 'color', color, 'LineWidth', 4);
        pd = weibull_pd;
        legendLabels = {'','Mean', 'Weibull fit', '2nd Quartile', '3rd Quartile'};
        %plot([ci_weibull(1), ci_weibull(1)], y_limits, '--k', 'LineWidth', 2, 'DisplayName', 'Lognormal 95% CI');
        %plot([ci_weibull(2), ci_weibull(2)], y_limits, '--k', 'LineWidth', 2);
        %text(ci_weibull(1), max(ylim) * 0.1, sprintf('%.2f', ci_weibull(1)), 'FontSize', 20, ...
        %    'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', 'Color', 'k');
        %text(ci_weibull(2), max(ylim) * 0.1, sprintf('%.2f', ci_weibull(2)), 'FontSize', 20, ...
        %    'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'Color', 'k');
        %legend({'','Mean','Median', 'Weibull fit', '95%CI', ''});
    else
        plot(x_values, pdf(lognormal_pd, x_values), 'color', color, 'LineWidth', 4);
        pd = lognormal_pd;
        legendLabels = {'','Mean', 'Lognormal fit', '2nd Quartile', '3rd Quartile'};
        %plot([ci_lognormal(1), ci_lognormal(1)], y_limits, '--k', 'LineWidth', 2, 'DisplayName', 'Lognormal 95% CI');
        %plot([ci_lognormal(2), ci_lognormal(2)], y_limits, '--k', 'LineWidth', 2);
        %text(ci_lognormal(1), max(ylim) * 0.1, sprintf('%.2f', ci_lognormal(1)), 'FontSize', 20, ...
        %    'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', 'Color', 'k');
        %text(ci_lognormal(2), max(ylim) * 0.1, sprintf('%.2f', ci_lognormal(2)), 'FontSize', 20, ...
        %    'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'Color', 'k');
        %legend({'','Mean','Median', 'Lognormal fit', '95%CI', ''});
    end

    ci_25 = icdf(pd, 0.25)
    ci_75 = icdf(pd, 0.75)

    % plot([ci_25, ci_25], y_limits, '--k', 'LineWidth', 2, 'DisplayName', 'Interquartile Range');
    % plot([ci_75, ci_75], y_limits, '--k', 'LineWidth', 2);

    % Shade the 2nd quartile (25th percentile to median) in light gray
    fill([ci_25, ci_25, medianC0, medianC0], [0, y_limits(end), y_limits(end), 0], ...
         [0.8, 0.8, 0.8], 'FaceAlpha', 0.4, 'EdgeColor', 'none');  % Light gray fill

    % Shade the 3rd quartile (median to 75th percentile) in slightly darker gray
    fill([medianC0, medianC0, ci_75, ci_75], [0, y_limits(end), y_limits(end), 0], ...
         [0.6, 0.6, 0.6], 'FaceAlpha', 0.4, 'EdgeColor', 'none');  % Darker gray fill

    text(ci_25, max(ylim) * 0.1, sprintf('%.2f', ci_25), 'FontSize', 20, ...
       'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', 'Color', 'k');
    text(ci_75, max(ylim) * 0.1, sprintf('%.2f', ci_75), 'FontSize', 20, ...
       'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'Color', 'k');
    legend(legendLabels);
    if ~isempty(plotTitle)
        title(plotTitle, 'FontSize', 36);
    end
end
