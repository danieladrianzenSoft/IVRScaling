function [C, CI_avg, CF_avg, CE_avg,CS_avg, Mass_dimless, N, Mass_remaining, mR_inner, CB, C_0, M_0, T, params] = solve_diffusion_5C(t, condition, drug, varargin)
p = inputParser;
p.addRequired('t', @(x) length(x)>=1);
p.addRequired('condition', @isstring);
p.addRequired('drug', @isstring);
p.addOptional('M_0', 0, @isscalar);
p.addOptional('C_0', 0, @isscalar);
p.addOptional('numSteps', 0, @isscalar);
p.addOptional('currentStep', 0, @isscalar);
p.addOptional('nonDimensional',false,@islogical)
p.parse(t,condition,drug,varargin{:});

isNonDimensional = p.Results.nonDimensional;
M_0 = p.Results.M_0;
C_0 = p.Results.C_0;
numSteps = p.Results.numSteps;
currentStep = p.Results.currentStep;

%%
%%%% VALIDATION FOR STEP ANALYSIS %%%%
if (contains(condition, 'step'))
    if numSteps == 0
        warning('step analysis selected but no numSteps provided. Setting numSteps=10 by default');
        numSteps = 10;
    end

    if currentStep == 0
        warning('currentStep not provided!')
        return;
    end
end

%%%%%% INITIALIZING PARAMETERS %%%%%%
%% Mass Transport Parameters

if strcmp(drug, 'hydrophobic')
    D_I = 1.3e-9;
    D_F = 5.7e-6;     % same value as palm gel paper
    D_E = 1.5e-7;     % from cornell study
    D_S = 1.5e-8;     % log lower than cornell study value
    phi_IF = .001;
    phi_FE = 1;
    phi_ES = .1;
    k_B = 1.19/(60^2); % Rate constant drug clearance in stroma 1 log higher than hydrophilic value because of lipid content in cell wall of blood vessels
    k_L = (.00866./(60^2))/2; % CHECK THIS, WE'RE USING EFDA VALUES HERE
elseif strcmp(drug, 'hydrophilic')
    D_I = 1.44e-9;
    D_F = 5.7e-6;    % from palm gel paper
    D_E = 6.1e-8;    % from palm gel paper
    D_S = 4.52e-7;   % from palm gel paper
    phi_IF = 1;
    phi_FE = 1;
    phi_ES = 1;
    k_B = .119/(60^2); % Rate constant drug clearance in stroma
    k_L = (.00866./(60^2))/2; % Rate constant blood clearance from halflife for EFDA
    %k_L = 141./(60^2); % from jing
end

if strcmp(condition, 'human_avg') == 1  %do this one
    h_F = .02;
    k_F = 1.22/(60^2);
    A_F = 100;
    h_E = .02;
    h_S = .15;
    Vb = 5000;
    R_outer = 2.7; % cm
    R_inner = 1.94; % cm
    Vs = A_F * h_S;
elseif strcmp(condition, 'human_large') == 1  %do this one
    h_F = .02;
    k_F = 1.22/(60^2);
    A_F = 120;
    h_E = .025;
    h_S = .2;
    Vb = 7000;
    R_outer = 2.7; % cm
    R_inner = 1.94; % cm
    Vs = A_F * h_S;
elseif strcmp(condition, 'human_small') == 1  %do this one
    h_F = .02;
    k_F = 1.22/(60^2);
    A_F = 80;
    h_E = .015;
    h_S = .1;
    Vb = 3000;
    R_outer = 2.7; % cm
    R_inner = 1.94; % cm
    Vs = A_F * h_S;
elseif (strcmp(condition, 'human_step_increasing_size') == 1)
    h_F = .02;
    k_F = 1.22/(60^2);
    A_F_vec = linspace(80,120,numSteps);
    A_F = A_F_vec(currentStep);
    h_E_vec = linspace(.015,.025,numSteps);
    h_E = h_E_vec(currentStep);
    h_S_vec = linspace(.1,.2,numSteps);
    h_S = h_S_vec(currentStep);
    Vb_vec = linspace(3000,7000,numSteps);
    Vb = Vb_vec(currentStep);
    R_outer = 2.7; % cm
    R_inner = 1.94; % cm
    Vs = A_F * h_S;
elseif (strcmp(condition, 'human_step_increasing_size_follicular') == 1)
    h_F = .01;
    k_F = 0.61/(60^2);
    A_F_vec = linspace(70,130,numSteps);
    A_F = A_F_vec(currentStep);
    h_E = .01;
    h_S_vec = linspace(.126,.154,numSteps);
    h_S = h_S_vec(currentStep);
    Vb_vec = linspace(3500,6500,numSteps);
    Vb = Vb_vec(currentStep);
    R_outer = 2.7; % cm
    R_inner = 1.94; % cm
    Vs = A_F * h_S;
elseif (strcmp(condition, 'human_step_increasing_size_midcycle') == 1)
    h_F = .02;
    k_F = 1.22/(60^2);
    A_F_vec = linspace(70,130,numSteps);
    A_F = A_F_vec(currentStep);
    h_E = .02;
    h_S_vec = linspace(.126,.154,numSteps);
    h_S = h_S_vec(currentStep);
    Vb_vec = linspace(3500,6500,numSteps);
    Vb = Vb_vec(currentStep);
    R_outer = 2.7; % cm
    R_inner = 1.94; % cm
    Vs = A_F * h_S;
elseif (strcmp(condition, 'human_step_increasing_size_luteal') == 1)
    h_F = .01;
    k_F = 0.61/(60^2);
    A_F_vec = linspace(70,130,numSteps);
    A_F = A_F_vec(currentStep);
    h_E = .01;
    h_S_vec = linspace(.126,.154,numSteps);
    h_S = h_S_vec(currentStep);
    Vb_vec = linspace(3500,6500,numSteps);
    Vb = Vb_vec(currentStep);
    R_outer = 2.7; % cm
    R_inner = 1.94; % cm
    Vs = A_F * h_S;
elseif strcmp(condition, 'human_avg_luteal') == 1  %do this one
    h_F = .01;
    k_F = 0;
    A_F = 100;
    h_E = .01;
    h_S = .15;
    Vb = 5000;
    R_outer = 2.7; % cm
    R_inner = 1.94; % cm
    Vs = A_F * h_S;
elseif strcmp(condition, 'human_large_luteal') == 1  %do this one
    h_F = .01;
    k_F = 0;
    A_F = 120;
    h_E = .0125;
    h_S = .2;
    Vb = 7000;
    R_outer = 2.7; % cm
    R_inner = 1.94; % cm
    Vs = A_F * h_S;
elseif strcmp(condition, 'human_small_luteal') == 1  %do this one
    h_F = .01;
    k_F = 0;
    A_F = 80;
    h_E = .01;
    h_S = .1;
    Vb = 3000;
    R_outer = 2.7; % cm
    R_inner = 1.94; % cm
    Vs = A_F * h_S;
elseif strcmp(condition, 'macaque_avg') == 1  %do this one
    h_F = .015;
    k_F = 1.22/(60^2);
    A_F = 40;
    h_E = .015;
    h_S = .1;
    Vb = 1000;
    R_outer = 1.25 ; % cm
    R_inner = .65; % cm
    Vs = A_F * h_S;
elseif strcmp(condition, 'macaque_large') == 1  %do this one
    h_F = .015;
    k_F = 1.22/(60^2);
    A_F = 48;
    h_E = .0188;
    h_S = .15;
    Vb = 1400;
    R_outer = 1.25 ; % cm
    R_inner = .65; % cm
    Vs = A_F * h_S;
elseif strcmp(condition, 'macaque_small') == 1  %do this one
    h_F = .015;
    k_F = 1.22/(60^2);
    A_F = 32;
    h_E = .0114;
    h_S = .075;
    Vb = 600;
    R_outer = 1.25 ; % cm
    R_inner = .65; % cm
    Vs = A_F * h_S;
elseif (strcmp(condition, 'macaque_step_increasing_size') == 1)
    h_F = .015;
    k_F = 1.22/(60^2);
    A_F_vec = linspace(32,48,numSteps);
    A_F = A_F_vec(currentStep);
    h_E_vec = linspace(.0114,.0188,numSteps);
    h_E = h_E_vec(currentStep);
    h_S_vec = linspace(.075,.15,numSteps);
    h_S = h_S_vec(currentStep);
    Vb_vec = linspace(600,1400,numSteps);
    Vb = Vb_vec(currentStep);
    R_outer = 1.25 ; % cm
    R_inner = .65; % cm
    Vs = A_F * h_S;
elseif (strcmp(condition, 'macaque_step_increasing_size_follicular') == 1)
    h_F = .01;
    k_F = 0.61/(60^2);
    A_F_vec = linspace(32,59,numSteps);
    A_F = A_F_vec(currentStep);
    h_E = .02;
    h_S_vec = linspace(.09,.11,numSteps);
    h_S = h_S_vec(currentStep);
    Vb_vec = linspace(700,1300,numSteps);
    Vb = Vb_vec(currentStep);
    R_outer = 1.25 ; % cm
    R_inner = .65; % cm
    Vs = A_F * h_S;
elseif (strcmp(condition, 'macaque_step_increasing_size_midcycle') == 1)
    h_F = .02;
    k_F = 1.22/(60^2);
    A_F_vec = linspace(32,59,numSteps);
    A_F = A_F_vec(currentStep);
    h_E = .03;
    h_S_vec = linspace(.09,.11,numSteps);
    h_S = h_S_vec(currentStep);
    Vb_vec = linspace(700,1300,numSteps);
    Vb = Vb_vec(currentStep);
    R_outer = 1.25 ; % cm
    R_inner = .65; % cm
    Vs = A_F * h_S;
elseif (strcmp(condition, 'macaque_step_increasing_size_luteal') == 1)
    h_F = .01;
    k_F = 0.61/(60^2);
    A_F_vec = linspace(32,59,numSteps);
    A_F = A_F_vec(currentStep);
    h_E = .02;
    h_S_vec = linspace(.9,.11,numSteps);
    h_S = h_S_vec(currentStep);
    Vb_vec = linspace(700,1300,numSteps);
    Vb = Vb_vec(currentStep);
    R_outer = 1.25 ; % cm
    R_inner = .65; % cm
    Vs = A_F * h_S;
elseif strcmp(condition, 'macaque_avg_luteal') == 1  %do this one
    h_F = .0075;
    k_F = 0;
    A_F = 40;
    h_E = .0075;
    h_S = .1;
    Vb = 1000;
    R_outer = 1.25 ; % cm
    R_inner = .65; % cm
    Vs = A_F * h_S;
elseif strcmp(condition, 'macaque_large_luteal') == 1  %do this one
    h_F = .0075;
    k_F = 0;
    A_F = 48;
    h_E = .0075;
    h_S = .15;
    Vb = 1400;
    R_outer = 1.25 ; % cm
    R_inner = .65; % cm
    Vs = A_F * h_S;
elseif strcmp(condition, 'macaque_small_luteal') == 1  %do this one
    h_F = .0075;
    k_F = 0;
    A_F = 32;
    h_E = .0057;
    h_S = .075;
    Vb = 600;
    R_outer = 1.25 ; % cm
    R_inner = .65; % cm
    Vs = A_F * h_S;
elseif (strcmp(condition, 'sheep_step_increasing_size_follicular') == 1)
    h_F = .01;
    k_F = 0.61/(60^2);
    A_F_vec = linspace(53,97,numSteps);
    A_F = A_F_vec(currentStep);
    h_E = .007;
    h_S_vec = linspace(.117,.143,numSteps);
    h_S = h_S_vec(currentStep);
    Vb_vec = linspace(2730,5070,numSteps);
    Vb = Vb_vec(currentStep);
    R_outer = 2.7; % cm
    R_inner = 1.94; % cm
    Vs = A_F * h_S;
elseif (strcmp(condition, 'sheep_step_increasing_size_midcycle') == 1)
    h_F = .02;
    k_F = 1.22/(60^2);
    A_F_vec = linspace(53,97,numSteps);
    A_F = A_F_vec(currentStep);
    h_E = .01;
    h_S_vec = linspace(.117,.143,numSteps);
    h_S = h_S_vec(currentStep);
    Vb_vec = linspace(2730,5070,numSteps);
    Vb = Vb_vec(currentStep);
    R_outer = 2.7; % cm
    R_inner = 1.94; % cm
    Vs = A_F * h_S;
elseif (strcmp(condition, 'sheep_step_increasing_size_luteal') == 1)
    h_F = .01;
    k_F = 0.61/(60^2);
    A_F_vec = linspace(53,97,numSteps);
    A_F = A_F_vec(currentStep);
    h_E = .007;
    h_S_vec = linspace(.117,.143,numSteps);
    h_S = h_S_vec(currentStep);
    Vb_vec = linspace(2730,5070,numSteps);
    Vb = Vb_vec(currentStep);
    R_outer = 2.7; % cm
    R_inner = 1.94; % cm
    Vs = A_F * h_S;
end

    %R_outer = 2.7; % cm
    %R_inner = 1.94; % cm

%% Geometric Parameters
% R_outer = 2.7; % cm
% R_inner = 1.94; % cm
mR_inner = (R_outer-R_inner)/2;
mR_outer = R_outer - mR_inner;
%SSA_star = 2./mR_inner_star;
%SSA = SSA_star * SSA_change;
%mR_inner = 2./(SSA);

W = 5; % cm
L = 10; % cm

A_IVR = (2*pi*mR_outer)*(2*pi*mR_inner); % surface area of IVR (cm^2)
if strcmp(condition, 'sheep luteal') == 1
    A_IVR = A_IVR/2;
end
V_IVR = pi*mR_inner^2*2*pi*mR_outer; % volume of IVR (cm^3)
%A_F = 100;   % area of fluid
h_I = mR_inner; % radius of IVR (cm)
%h_S = .3; % thickness of stroma layer (cm)
%Vb = 7500;   %volume of blood
heights = [h_I, h_F+h_I, h_E+h_I+h_F, h_S+h_E+h_F+h_I];
%% Solver Parameters
numX = 3001; % total number of spatial points.
x = linspace(0,h_I+h_F+h_E+h_S,numX); % [cm] spatial vector, set to the total length of domain (IVR + tissue).
dx = x(2)-x(1);
nI = ceil(((h_I)/(h_I+h_F+h_E+h_S))*numX); %number of spatial points in the IVR
nF = ceil(((h_F)/(h_I+h_F+h_E+h_S))*numX); %number of spatial points in the epithelium
nE = ceil(((h_E)/(h_I+h_F+h_E+h_S))*numX);
nS = numX - (nI+nF+nE); %number of spatial points in stroma

if (isNonDimensional == true)
    if ((C_0 > 0 && C_0 ~= 1) || M_0 > 0 )
        warning('NonDimensional is set to true, but C_0 or M_0 inputs were provided. These inputs will be ignored, C_0 is set to 1.');
    end
    C_0 = 1;
    M_0 = C_0*V_IVR;
elseif (C_0 > 0)
    if (M_0 > 0)
       warning('Both M_0 and C_0 were provided. C_0 is set to %.2f and M_0 input is ignored', C_0);
    end
    M_0 = C_0*V_IVR;
elseif (M_0 > 0)
    C_0 = M_0/V_IVR;
elseif (C_0 <= 0 && M_0 <= 0)
    warning('M_0 and C_0 are invalid at %.2f and %.2f respectively', M_0, C_0);
end

%EC50 = 50.*287; %ng/mL
EC50 = 2.93e-2;
%.1 nM -> 2.93e-8 mg/mL
% 1 nM -> 2.93e-7 mg/mL
% 5 nM -> 1.47e-6 mg/mL
%% Initial Condition
IC = zeros(1,nI+nF+nE+nS);
IC(1:nI) = C_0;

%% Run the Model
t_og = t;
tic
[t,C] = ode15s(@(t2,C) Diffusion_1D_4Cmpts(t2,C,x,D_I,D_F,D_E,D_S,k_F,k_B, phi_IF,phi_FE,phi_ES, nI,nF,nE,nS, A_F, A_IVR), t, IC);
toc
C = C';
t_hr = t./3600; %convert to hr

%Cbiopsy_1 = C(nI+nF+1:biopsyIndex, :);
%Cf_end_1 = Cf(:,end);
%Ce_end_1 = Ce(:,end);

%% Reduce SA
% [t,C_red] = ode15s(@(t2,C_red) Diffusion_1D_4Cmpts(t2,C_red,x,D_I,D_F,D_E,D_S,k_F,k_B, phi_IF,phi_FE,phi_ES, nI,nF,nE,nS, A_F, A_IVR./2), t, IC);
% C_red = C_red';
% 
% Ci_red = C_red(1:nI,:);
% Cf_red = C_red(nI+1:nI+nF,:);
% Ce_red = C_red(nI+nF+1:nI+nF+nE,:);
% Cs_red = C_red(nI+nF+nE+1:end, :);

%% Lag IC
% IC_lag = C(nI+1:end, end)';
% t3 = linspace(5*24*60*60+1, 28*24*60*60,10000);
% x2 = linspace(0, h_F+h_E+h_S,numX-nI);
% [t3,C2] = ode15s(@(t4,C2) Diffusion_1D_4Cmpts_Lag(t3,C2,x2,D_I,D_F,D_E,D_S,k_F,k_B, phi_IF,phi_FE,phi_ES, nI,nF,nE,nS, A_F, A_IVR), t3, IC_lag);
% t = [t; t3];
% C2 = C2';
% Cf_2 = C2(1:nF,:);
% Ce_2 = C2(nF+1:nF+nE,:);
% Cs_2 = C2(nF+nE+1:end, :);
% 
% Cf = horzcat(Cf_1, Cf_2);
% Ce = horzcat(Ce_1, Ce_2);
% Cs = horzcat(Cs_1, Cs_2);

%% Simulated Biopsy
% Biopsy depth = 2mm
%biopsyIndex = find(x>=.2+h_F+h_I, 1, 'first');
%nBiopsy = biopsyIndex - (nI+nF);
%% Average Concentrations
Ci = C(1:nI,:);
Cf = C(nI+1:nI+nF,:);
Ce = C(nI+nF+1:nI+nF+nE,:);
Cs = C(nI+nF+nE+1:end, :);
%Cbiopsy = C(nI+nF+1:biopsyIndex, :);
Cf_end = Cf(:,end);
Ce_end = Ce(:,end);
name = strcat(condition,'.mat');
%save(name,"x", "Cf_end", "Ce_end");

%CBiopsy_avg = trapz(Cbiopsy)./(nBiopsy-1);
CI_avg = trapz(Ci)/(nI-1);
CF_avg = trapz(Cf)/(nF-1);
CE_avg = trapz(Ce)/(nE-1);
CS_avg = trapz(Cs)/(nS-1);
C_nonring = C(nI+1:end, :);

% CI_avg_red = trapz(Ci_red)/(nI-1);
% CF_avg_red = trapz(Cf_red)/(nF-1);
% CE_avg_red = trapz(Ce_red)/(nE-1);
% CS_avg_red = trapz(Cs_red)/(nS-1);



%% Flux and Mass Release
%Mass = trapz(x(nI+1:end), C_nonring) .* W .* L;
%N = D_I.* abs(C(nI+1,:)-C(nI,:))./dx;
N = -1.* D_I.*(C(nI,:) - C(nI-1,:)) ./ dx;  % compute flux at the edge of the slab
Mass_dot = N .* A_IVR;                           % compute mass released as the integral of flux
Mass = cumtrapz(t_og, Mass_dot);
Mass_dimless = Mass./M_0;

%% Residual Drug Computation
new_C = zeros(size(Ci));
for i = 1:nI
    for j = 1:length(t_og)
        new_C(i,j) = Ci(i,j)*2*pi*x(i);
    end
end
Mass_remaining = (2*pi*mR_outer*dx*trapz(new_C)+ dx*A_F*trapz(C(nI+1:end, :)))./M_0;

if (0)
    figure;
    plot(t_og, Mass_remaining)
end

%% Blood
CB = zeros(1,length(t));
%f = CS_avg*W*L*h_S*k_B*2;
%decay = exp(-t*k_L)';
lossinblood = zeros(1, length(t_hr));
V_S = h_S.*W*L;
dt = t(2)-t(1);
for i = 2:length(t)
    dCBdT = ((CS_avg(i-1).*k_B.*V_S - CB(i-1).*k_L.*Vb)./Vb).*dt;
    CB(i) = CB(i-1)+dCBdT;
    %CB(i) = trapz(t(1:i),f(1:i).*decay(i:-1:1)/Vb);
    lossinblood(i) = CB(i-1).*k_L;
end

%% Blood for reduced
% CB_red = zeros(1,length(t));
% V_S = h_S.*W*L;
% dt = t(2)-t(1);
% for i = 2:length(t)
%     dCBdT = ((CS_avg_red(i-1).*k_B.*V_S - CB_red(i-1).*k_L.*Vb)./Vb).*dt;
%     CB_red(i) = CB_red(i-1)+dCBdT;
% end
%% Mass in each Compartment

lossinfluid = CF_avg.*k_F;

%MassRing = CI_avg.*V_IVR;
%MassFluid = CF_avg.*h_F.*A_F;
%MassEpithelium = CE_avg.*h_E.*A_F;
%MassStroma = CS_avg.*h_S.*A_F;
MassRing = 2*pi*mR_outer*dx*trapz(new_C);
MassFluid = dx.*A_F.*trapz(C(nI+1:nI+nF, :));
MassEpithelium = dx.*A_F.*trapz(C(nI+nF+1:nI+nF+nE, :));
MassStroma = dx.*A_F.*trapz(C(nI+nF+nE+1:nI+nF+nE+nS, :));
MassBlood = CB.*Vb;
MassLossFluid = cumtrapz(t, lossinfluid).*h_F.*A_F;
MassLossBlood = cumtrapz(t, lossinblood).*Vb;
TotalMass = MassRing + MassFluid + MassEpithelium + MassStroma + MassBlood + MassLossFluid + MassLossBlood;
%% Animated Pie Chart
if 0
    % Each row is a time point, each column is a category
    data = [MassRing; MassFluid; MassEpithelium; MassStroma; MassBlood; MassLossFluid; MassLossBlood];
    %data = data ./ sum(data, 2); % Normalize to ensure each row sums to 1
    data = data(:,1:100:end);
    t_sampled = t(1:100:end);
    % Number of time points
    numTimePoints = size(data, 2);
    labels = {'Ring', 'Fluid', 'Epithelium', 'Stroma', 'Blood','Fluid Loss', 'Blood Loss'};
    % Prepare the figure
    figure;

    % Create a VideoWriter object
    videoFile = strcat(condition, '.mp4');
    v = VideoWriter(videoFile, 'MPEG-4');
    v.FrameRate = 20; % Set frame rate
    open(v);


    % Loop over each time point
    for  p = 1:numTimePoints
        % Draw the pie chart for the current time point
        clf; % Clear the current figure
        pie(data(:, p));
        txt= sprintf('%.1f Days', t_sampled(p)./(60*60*24));
        text(-1, 1, txt, 'FontSize', 12); % Position the text
        legend(labels, 'Location', 'bestoutside');

        % Capture the frame
        frames(p) = getframe(gcf);

        % Write frame to video
        writeVideo(v, frames(p));

    end
    close(v);

    % Play back the animation
    %movie(frames, 1, 100); % Play once at 10 frames per second
end
%% Plots of Mass in each compartment
if 0
    figure;
    plot(t./(60*60*24), MassRing, '-k', 'LineWidth', 2, 'DisplayName', 'Ring')
    hold on
    plot(t./(60*60*24), MassFluid, '-', 'Color', [0.3010 0.7450 0.9330], 'LineWidth', 2, 'DisplayName', 'Fluid')
    plot(t./(60*60*24), MassEpithelium, '-', 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 2, 'DisplayName', 'Epithelium')
    plot(t./(60*60*24), MassStroma, '-', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2, 'DisplayName', 'Stroma')
    plot(t./(60*60*24), MassBlood, '-r', 'LineWidth', 2, 'DisplayName', 'Blood')
    plot(t./(60*60*24), MassLossFluid, '-', 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 2, 'DisplayName', 'Fluid Loss')
    plot(t./(60*60*24), MassLossBlood, '-', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2, 'DisplayName', 'Blood Loss')
    xlabel('time [day]', 'FontSize', 15)
    ylabel('Mass in Compartment', 'FontSize',15)
    title(strcat(condition, ' Mass Partitioning'), 'FontSize', 20)
    set(gca, 'FontSize', 20)
    legend show

    figure;
    plot(t./(60*60*24), TotalMass, '--k', 'LineWidth',2)
    xlabel('time [days]', 'FontSize',15)
    ylabel(strcat('Mass in System',condition), 'FontSize',15)
    set(gca, 'FontSize', 20)


end
idx_1 = find(t>=1*24*60*60, 1);
idx_3 = find(t>=3*24*60*60, 1);
idx_8 = find(t>=8*24*60*60, 1);
Plasma_3_1 = CB(idx_3)./CB(idx_8);
Plasma_8_3 = CB(idx_8)./CB(idx_3);

Fluid_3_1 = CF_avg(idx_3)./CF_avg(idx_1);
Fluid_8_3 = CF_avg(idx_8)./CF_avg(idx_3);


%% PP
if 0
    prot_matrix = C>EC50;
    S = sum(prot_matrix(nI+nF+1:end,:),1);

    PP = S./(nE+nS);
    plateau = diff(PP);

    tlag = t_hr(find(PP == max(PP), 1, 'first'))
    PPmax = PP(find(PP == max(PP), 1, 'first'))

    figure(10);
    plot(t./(60*60*24),PP,'LineWidth', 2, 'DisplayName', 'Dimensionless')
    hold on
    ylabel('Percent of Tissue Protected (%)', 'FontSize',15)
    xlabel('time [day]', 'FontSize',15)
    title('PP', 'FontSize',20)
    legend show
end

%save('concentration_smallC0.mat', "t", "CI_avg", "CF_avg", "CE_avg", "CS_avg", "CB", "CBiopsy_avg", "Mass_remaining")

%% Plots

if 0
    AvgConcStacked(t, CI_avg, CF_avg, CE_avg, CS_avg, D_I, condition, drug, CB)    
    figure;
    plot(t./(60*60*24), Mass_remaining, '-k', 'LineWidth', 2)
    xlabel('time [day]', 'FontSize',15)
    ylabel('Mass in Ring [mg]', 'FontSize',15)
    
    figure;
    plot(t./(60*60*24), Mass_dimless, '-k', 'LineWidth', 2)
    xlabel('time [day]', 'FontSize',15)
    ylabel('Mass Released [mg]', 'FontSize',15)
end

if 0
    Depth(C, x, t, condition, drug, heights, 4, nI)
end

%% PK Data
strcat(drug, ', ', condition);
params = table(h_F, k_F, h_E, h_S, Vb, R_outer, R_inner, Vs);

T = ParAnalysis(t, CI_avg, CF_avg, CE_avg, CS_avg, Mass_dimless, CB, C_0);
T = rows2vars(T,"VariableNamingRule","preserve");

end

