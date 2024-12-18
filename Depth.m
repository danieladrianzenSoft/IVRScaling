function Depth(C,x,t, cyclePhase, drug, heights, num_compartments, nI)

 titlestr = strcat(drug,{' '}, 'in', {' '}, cyclePhase);

figure;
hours =  [1 12 24 24*7 24*14 ];

color1 = [0.8 0.8 1]; % light blue
color2 = [1 0.8 0.8]; % light red
color3 = [1 0.9 0.6]; % light orange
color4 = [0.9 0.8 0.8]; % light pink
patch([0 heights(1) heights(1) 0], [0 0 max(ylim) max(ylim)], color1, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
patch([heights(1) heights(2) heights(2) heights(1)], [0 0 max(ylim) max(ylim)], color2, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
patch([heights(2) heights(3) heights(3) heights(2)], [0 0 max(ylim) max(ylim)], color3, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

if num_compartments == 4
patch([heights(3) heights(4) heights(4) heights(3)], [0 0 max(ylim) max(ylim)], color4, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
end

for i =1:length(hours)
    hour_names{i} = sprintf('%d hours', hours(i));
end
for i=hours
    plot(x(1:end)',C(1:end,find(t >= i*60^2,1,'first')), 'LineWidth', 2)
    hold on
end

for i = 1:length(heights)
    xline(heights(i), '--k')
end
%xlim([ .38 .7])
hold off
set(gca,'FontSize',20)
legend(hour_names,'Location','Northeast','FontSize',20)
ylabel("Concentration (mg/mL)",'FontSize',24)
xlabel("Depth (cm)",'FontSize',24)
title(titlestr,'FontSize',24)
filename = strcat('concprof_', drug, '_', cyclePhase, '.png');
saveas(gcf, filename)

end