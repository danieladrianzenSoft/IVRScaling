function AvgConcStacked(t,CI_avg, CF_avg, CE_avg,CS_avg, D_I, cyclePhase, drug, CB)
td = t/(60^2*24);
% titlestr = strcat('Averge Concentration for ',{' '} ,loading,{' '}, 'in', {' '}, cyclePhase);

legend_str = strcat(cyclePhase);

if strcmp(cyclePhase, 'sheep luteal') == 1 
linspec = '-k';
elseif strcmp(cyclePhase, 'sheep midcycle') == 1  %do this one
linspec = '-r';
elseif strcmp(cyclePhase, 'macaque luteal') == 1  %do this one
linspec = '-k';
elseif strcmp(cyclePhase, 'macaque midcycle') == 1  %do this one
linspec = '-r';
else
    linspec = '-b';

end

figure(55)
subplot(2, 2, 1)
plot(td,CE_avg,linspec,'LineWidth', 2,'DisplayName', cyclePhase)
legend show
hold on
%ylabel("Concentration",'FontSize',15)
xlabel("Time Post Insertion (day)",'FontSize',15)
ylabel('Concentration [mg/mL]', 'FontSize',15)
title('Vaginal Epithelium', 'FontSize',17)
set(gca, 'FontSize', 17)

subplot(2, 2, 2)
hold on
plot(td,CS_avg,linspec,'LineWidth',2,'DisplayName', legend_str)
hold on
%ylabel("Concentration",'FontSize',15)
xlabel("Time Post Insertion (day)",'FontSize',15)
ylabel('Concentration [mg/mL]', 'FontSize',15)
title('Vaginal Stroma', 'FontSize',17)
set(gca, 'FontSize', 17)


subplot(2, 2, 3)
hold on
plot(td, CF_avg, linspec,'LineWidth', 2,'DisplayName', legend_str)
hold on
%ylabel("Concentration",'FontSize',15)
xlabel("Time Post Insertion (day)",'FontSize',15)
ylabel('Concentration [mg/mL]', 'FontSize',15)
title({'';'Vaginal Luminal Fluid'}, 'FontSize',17)
set(gca, 'FontSize', 17)


subplot(2, 2, 4)
hold on
plot(td, CB,linspec, 'LineWidth',2,'DisplayName',legend_str)
hold on
%ylabel("Concentration",'FontSize',15)
xlabel("Time Post Insertion (day)",'FontSize',17)
ylabel('Concentration [mg/mL]', 'FontSize',15)
title({'';'Blood'}, 'FontSize',17)
set(gca, 'FontSize', 17)

% subplot(5, 1, 5)
% hold on
% plot(td, CBiopsy_avg, 'LineWidth', 2,'DisplayName', legend_str)
% hold on
% legend show
% ylabel("Concentration",'FontSize',15)
% xlabel("Time (day)",'FontSize',15)
% title('Simulated Biopsy Concentration Plots', 'FontSize',17)
% xlim([0 30])

end