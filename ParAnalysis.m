function [T] = ParAnalysis(t,CI_avg,CF_avg, CE_avg,CS_avg, mass_dimless, CB, C_0)

Compt = {CF_avg CF_avg CE_avg CS_avg CB};
n = length(Compt);
slope1 = cell(1,n);
slope2 = cell(1,n);
abslope2 = cell(1,n);
imax = zeros(1,n);
platI = zeros(1,n);
Cmax = zeros(1,n);
Tmax = zeros(1,n);
C14 = zeros(1,n);
C28 = zeros(1,n);
C90 = zeros(1,n);
AUC14 = zeros(1,n);
AUC28 = zeros(1,n);
AUC90 = zeros(1,n);
Cstar28 = zeros(1,n);
itmax = zeros(1,n);
AUCtmax = zeros(1,n);

Tplat = zeros(1,n);
Cplat = zeros(1,n);
R = zeros(1,n);

%i1 = find(t >= 1*60^2*24,1,'first');
%i2 = find(t >= 2*60^2*24,1,'first');
%i7 = find(t >= 7*60^2*24,1,'first');
i14 = find(t >= 14*60^2*24,1,'first');
i28 = find(t >= 28*60^2*24,1,'first');
%i90 = find(t >= 90*60^2*24,1,'first');
Mass_14 = mass_dimless(i14);
%Mass_28 = mass_dimless(i28);
%Mass_90 = mass_dimless(i90);

blankArray = zeros(1, n-1);
Mass = [Mass_14];
Mass = [Mass, blankArray];
for i = 1:length(Compt)
    slope1{i} = gradient(Compt{i})'./gradient(t);
    slope2{i} = gradient(slope1{i})./gradient(t);
    abslope2{i} = abs(slope2{i});
    imax(i) = find(Compt{i} == max(Compt{i}));
    Cmax(i) = Compt{i}(imax(i));
    Tmax(i) = t(imax(i))./60^2;
    itmax(i) = find(t >=Tmax.*(60^2), 1, 'first');
%     C28(i) = Compt{i}(i28);
%     AUC28(i) = trapz(t./60^2,Compt{i});
%     C1(i) = Compt{i}(i1);
%     AUC1(i) = trapz(t(1:i1)./60^2,Compt{i}(1:i1));
%     C2(i) = Compt{i}(i2);
%     AUC2(i) = trapz(t(1:i2)./60^2,Compt{i}(1:i2));
%     C7(i) = Compt{i}(i7);
%     AUC7(i) = trapz(t(1:i7)./60^2,Compt{i}(1:i7));
    C14(i) = Compt{i}(i14);
    AUC14(i) = trapz(t(1:i14)./60^2,Compt{i}(1:i14));
    C28(i) = Compt{i}(i28);
   % C90(i) = Compt{i}(i90);
    AUC28(i) = trapz(t(1:i28)./60^2,Compt{i}(1:i28));
    Cstar28(i) = AUC28(i)/(28*24);
   % AUC90(i) = trapz(t(1:i90)./60^2,Compt{i}(1:i90));
   % AUCtmax(i) = trapz(t(1:imax(i))./60^2,Compt{i}(1:imax(i)));
%     [M,platI(i)] = min(abs(Compt{i}(i28)./Compt{i}-0.2));
%     Tplat(i) = t(platI(i))./60^2;
%     Cplat(i) = Compt{i}(platI(i));
%     R(i) = 2*(Cplat(i)-C28(i))./(Cplat(i)+C28(i));
end

%par = [Cmax; Tmax; C28; AUC28; Tplat; Cplat; R];
%par = [Cmax; Tmax; C28; AUC28; C1; AUC1; C2; AUC2; C7; AUC7];

%par = [AUCtmax; Cmax; Tmax; C28; AUC28; Mass];
par = [Tmax;Cmax/C_0; C14/C_0; AUC14/C_0; C28/C_0; AUC28/C_0; Cstar28/C_0; Mass];
% T = array2table(par, 'VariableNames',...
%     {'IVR','Epithelium', 'Stroma', 'Blood'}, 'RowNames',...
%     {'Cmax (ng/mL)', 'Tmax (hr)', 'C28 (ng/mL)', 'AUC28 (ng*hr/mL)','Tpl (hr)', 'Cpl (ng/mL)', 'R (a.u.)'});
% T = array2table(par, 'VariableNames',...
%     {'IVR','Epithelium', 'Stroma', 'Blood'}, 'RowNames',...
%     {'Cmax (ng/mL)', 'Tmax (hr)', 'C28 (ng/mL)', 'AUC28 (ng*hr/mL)',...
%     'C1 (ng/mL)', 'AUC1 (ng*hr/mL)', 'C2 (ng/mL)', 'AUC2 (ng*hr/mL)',...
%     'C7 (ng/mL)', 'AUC7 (ng*hr/mL)'});

% T = array2table(par, 'VariableNames',...
%     {'IVR','Fluid', 'Epithelium', 'Stroma', 'Blood', 'Simulated Biopsy'}, 'RowNames',...
%     {'Tmax (hr)', 'Cmax (mg/mL)', 'C14 (mg/mL) ', 'AUC14', 'Mass Release'});

T = array2table(par, 'VariableNames',...
    {'IVR','Fluid', 'Epithelium', 'Stroma', 'Blood'}, 'RowNames',...
    {'Tmax (hr)', 'Cmax/C0', 'C14/C0', 'AUC14/C0 (hrs)', 'C28/C0', 'AUC28/C0 (hrs)', 'C*28/C0', 'Mass Release'});


end