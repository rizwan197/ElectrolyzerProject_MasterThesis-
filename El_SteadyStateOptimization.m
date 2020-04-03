function [z0, x0, u0,EXIT] = El_SteadyStateOptimization(N,X0,P0)

import casadi.*
par = parElectrolyzer(N);

%% Build the plant model and solve steady state optimization problem
[xDiff, xAlg, input, eqnAlg, eqnDiff] = model(par.N);
x = [xAlg;xDiff]; 

%% Defining the disturbance
Pnet = SX.sym('Pnet');
Ptot = SX.zeros(1,1);

for nEl = 1:par.N
    Ptot = Ptot + xAlg(2*par.N+nEl);
end
eqnPnet = Pnet - Ptot;%Total power = sum of power of the individual electrolyzer

%% preparing symbolic variables
w = {};
% preparing numeric variables and bounds
w0 = [];
lbw = [];
ubw = [];

% declaring them symbolic
w = {w{:},x,input};
% declaring them numerically
w0 = [w0;X0']; %initial guess

%defining constraints on the decision variables (states and inputs i.e. MVs)

%constraints on states
lbu_k = 0*ones(par.N,1);%lower bound on cell voltage
ubu_k = inf*ones(par.N,1);

lbi_k = zeros(par.N,1);%lower bound on current
ubi_k = inf*ones(par.N,1);

lbP_k = 0*ones(par.N,1);%lower bound on power
ubP_k = inf*ones(par.N,1);

lbFeff_k = 0*ones(par.N,1);%lower bound on faraday efficiency
ubFeff_k = 1*ones(par.N,1);

lbnH2_k = zeros(par.N,1);%lower bound on hydrogen production
ubnH2_k = inf*ones(par.N,1);

lbqH2Oloss_k = zeros(par.N,1);%lower bound on water loss
ubqH2Oloss_k = inf*ones(par.N,1);

lbnH2el_net = 0;%lower bound on net hydrogen production from the electrolyzer
ubnH2el_net = inf;

lbnH2out_net = 0;%lower bound on hydrogen from the storage tank outlet 
ubnH2out_net = inf;

lbnO2el_net = 0;%lower bound on net oxygen production from the electrolyzer
ubnO2el_net = inf;

lbnO2out_net = 0;%lower bound on oxygen from the storage tank outlet
ubnO2out_net = inf;

lbT_el_out = 0;%lower bound on the temperature at electrolyzer outlet
ubT_el_out = inf;

lbT_k = 25*ones(par.N,1);%lower bound on the electrolyzer temperature
ubT_k = 80*ones(par.N,1);

lbPstoH2 = 20;%lower bound on the hydrogen storage pressure
ubPstoH2 = 30;

lbPstoO2 = 20;%lower bound on the oxygen storage pressure
ubPstoO2 = 30;

lbMbt = 0;%lower bound on the mass in the buffer tank 
ubMbt = 6000000;

lbT_bt_out = 0;%lower bound on the temperature of lye leaving the buffer tank 
ubT_bt_out = inf;

lbT_el_in = 0;%lower bound on the temperature at electrolyzer inlet
ubT_el_in = inf;

lbT_cw_out = 0;%lower bound on the coolant outlet temperature
ubT_cw_out = inf;

%constraints on the inputs
lbU_el_k = zeros(par.N,1);%lower bound on the electrolyzer voltage
ubU_el_k = inf*ones(par.N,1);
lbq_lye_k = 500*ones(par.N,1);%lower bound on the lye flowrate
ubq_lye_k = 10000*ones(par.N,1);
lbq_cw = 0.01;%lower bound on the coolant flow rate
ubq_cw = 80000;
lbzH2 = 0;%lower bound on hydrogen outlet valve opening
ubzH2 = 1;
lbzO2 = 0;%lower bound on oxygen outlet valve opening
ubzO2 = 1;
lbqH2O = 0;%lower bound on total water lost during electrolysis
ubqH2O = inf;

lbw = [lbw;lbu_k;lbi_k;lbP_k;lbFeff_k;lbnH2_k;lbqH2Oloss_k;lbnH2el_net;...
    lbnH2out_net;lbnO2el_net;lbnO2out_net;lbT_el_out;lbT_k;lbPstoH2;...
    lbPstoO2;lbMbt;lbT_bt_out;lbT_el_in;lbT_cw_out;lbU_el_k;lbq_lye_k;...
    lbq_cw;lbzH2;lbzO2;lbqH2O];%bounds on all the variables
ubw = [ubw;ubu_k;ubi_k;ubP_k;ubFeff_k;ubnH2_k;ubqH2Oloss_k;ubnH2el_net;...
    ubnH2out_net;ubnO2el_net;ubnO2out_net;ubT_el_out;ubT_k;ubPstoH2;...
    ubPstoO2;ubMbt;ubT_bt_out;ubT_el_in;ubT_cw_out;ubU_el_k;ubq_lye_k;...
    ubq_cw;ubzH2;ubzO2;ubqH2O]; 
 

%% preparing symbolic constraints
g = {};
% preparing numeric bounds
lbg = [];
ubg = [];

% declaring constraints
uElconst = [];
for nEl=1:par.N-1
    uElconst = [uElconst;xAlg(nEl)-xAlg(nEl+1)];
end

deltaT1 = xDiff(par.N+5)-par.Tw_in;
deltaT2 = xDiff(par.N+4)-xDiff(par.N+6);

deltaT_El1 = xDiff(1) - xDiff(par.N+5);%impose the connstraint on the maximum temperature difference between the inlet lye temp. and the electrolyzer temp.
deltaT_El2 = xDiff(2) - xDiff(par.N+5);
deltaT_El3 = xDiff(3) - xDiff(par.N+5);

Iden = SX.zeros(par.N,1);
for nEl = 1:par.N
    Iden(nEl) = (0.1*xAlg(par.N+nEl))/par.EL(nEl).A; %current density in mA/cm2
end
IdenMin = 32;   %minimum current density, 32 mA/cm2
IdenMax = 198.5;%maximum current density, 198.5 mA/cm2

g = {g{:},eqnAlg, eqnDiff,uElconst, Iden,eqnPnet,deltaT_El1,deltaT_El2,deltaT_El3,deltaT1,deltaT2};
lbg = [lbg;zeros(7*par.N+11,1);zeros(par.N-1,1);IdenMin*ones(par.N,1);0;0;0;0;2e-3;2e-3];
ubg = [ubg;zeros(7*par.N+11,1);zeros(par.N-1,1);IdenMax*ones(par.N,1);P0;30;30;30;inf;inf];

% g = {g{:},eqnAlg, eqnDiff,uElconst, Iden,eqnPnet};
% lbg = [lbg;zeros(7*par.N+11,1);zeros(par.N-1,1);IdenMin*ones(par.N,1);0];
% ubg = [ubg;zeros(7*par.N+11,1);zeros(par.N-1,1);IdenMax*ones(par.N,1);P0];

Objvol_H2 = SX.zeros(par.N,1);
for nEl = 1:par.N
    Objvol_H2(nEl) = (xAlg(4*par.N+nEl)*0.0224136*3600);%[Nm3/h]
end
J = -(Objvol_H2(1)+Objvol_H2(2)+Objvol_H2(3));

% J = 10;

%% formalize into an NLP problem
nlp = struct('x',vertcat(w{:}),'g',vertcat(g{:}),'f',J,'p',Pnet);

options = struct;
options.ipopt.print_level = 0;

% assign solver - IPOPT in this case
solver = nlpsol('solver','ipopt',nlp,options);

% solve - using the defined initial guess and bounds
sol = solver('x0',w0,'lbx',lbw,'ubx',ubw,'lbg',lbg,'ubg',ubg,'p',P0);
res = full(sol.x);

EXIT = solver.stats.return_status; 

%% Extracting results
Uk = [];
Ik = [];
Pk = [];
Feffk = [];
nH2k = [];
qH2Olossk = [];
Tk = [];
Vss = [];
q_lyek = [];

for nEl = 1:par.N
    %optimal value of the algebriac state
    Uk = [Uk res(nEl)];                         %cell voltage of the electrolyzer
    Ik = [Ik res(par.N+nEl)];                   %current in the electrolyzer
    Pk = [Pk res(2*par.N+nEl)];                 %power of the individual electrolyzer
    Feffk = [Feffk res(3*par.N+nEl)];           %faraday efficiency of each electrolyzer
    nH2k = [nH2k res(4*par.N+nEl)];             %hydrogen produced form each individual electrolyzer
    qH2Olossk = [qH2Olossk res(5*par.N+nEl)];   %water loss during electrolysis in kth electrolyzer
    %optimal value of the differential state
    Tk = [Tk res(6*par.N+5+nEl)];               %temperature of the individual electrolyzer
    %optimal value of the inputs
    Vss = [Vss res(7*par.N+11+nEl)];            %electrolyzer voltage
    q_lyek = [q_lyek res(8*par.N+11+nEl)];      %lye flowrate
end

%optimal value of the algebriac state
nH2El_tot=res(6*par.N+1);
nH2out_tot=res(6*par.N+2);
nO2El_tot=res(6*par.N+3);
nO2out_tot=res(6*par.N+4);
T_el_out=res(6*par.N+5);        %temp after mixing of the exiting liquid streams from all electrolyzers
%optimal value of the differential state
PstoH2=res(7*par.N+6);
PstoO2=res(7*par.N+7);
massBt=res(7*par.N+8);
T_bt_out=res(7*par.N+9);
T_el_in=res(7*par.N+10);         %temp of inlet lye stream coming into the electrolyzer
T_CW_out=res(7*par.N+11);
%optimal value of the inputs
qf_cw=res(9*par.N+12);
zH2=res(9*par.N+13);
zO2=res(9*par.N+14);
Qwater=res(9*par.N+15);


%% Calculation of initial state vector

%Nominal load H2 production and specific electricity consumption
V_H2_ini = nH2k*0.0224136*3600;%[Nm3/h]
for nEl = 1:par.N
    Ps_ini(nEl) = (Uk(nEl)*Ik(nEl)*par.EL(nEl).nc)/(1000*V_H2_ini(nEl));%[kWh/Nm3]
    Iden(nEl) = 0.1*Ik(nEl)/par.EL(nEl).A; 
end

Pnet = sum(Pk);
Ps_ini;
Eff_El = 3.55./Ps_ini;

z0 = [Uk Ik Pk Feffk nH2k qH2Olossk nH2El_tot nH2out_tot nO2El_tot nO2out_tot T_el_out];
x0 = [Tk PstoH2 PstoO2 massBt T_bt_out T_el_in T_CW_out];
u0 = [Vss q_lyek qf_cw zH2 zO2 Qwater];

end