function [z0, x0, u0] = El_SteadyStateOptimization(N,X0,P0)

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

lbT_k = 25*ones(par.N,1);%lower bound on the electrolyzer temperature
ubT_k = 80*ones(par.N,1);

lbMbt = 0*ones(par.N,1);%lower bound on the mass in the buffer tank 
ubMbt = 2000000*ones(par.N,1);

lbT_bt_out = 0*ones(par.N,1);%lower bound on the temperature of lye leaving the buffer tank 
ubT_bt_out = inf*ones(par.N,1);

lbT_el_in = 0*ones(par.N,1);%lower bound on the temperature at electrolyzer inlet
ubT_el_in = inf*ones(par.N,1);

lbT_cw_out = 0*ones(par.N,1);%lower bound on the coolant outlet temperature
ubT_cw_out = inf*ones(par.N,1);

%constraints on the inputs
lbU_el_k = zeros(par.N,1);%lower bound on the electrolyzer voltage
ubU_el_k = inf*ones(par.N,1);
lbq_lye_k = 500*ones(par.N,1);%lower bound on the lye flowrate
ubq_lye_k = 8000*ones(par.N,1);
lbq_cw = 1e-2*ones(par.N,1);%lower bound on the coolant flow rate
ubq_cw = 20000*ones(par.N,1);
lbqH2O = 0*ones(par.N,1);%lower bound on total water lost during electrolysis
ubqH2O = inf*ones(par.N,1);

lbw = [lbw;lbu_k;lbi_k;lbP_k;lbFeff_k;lbnH2_k;lbqH2Oloss_k;lbT_k;lbMbt;...
    lbT_bt_out;lbT_el_in;lbT_cw_out;lbU_el_k;lbq_lye_k;lbq_cw;lbqH2O];%bounds on all the variables
ubw = [ubw;ubu_k;ubi_k;ubP_k;ubFeff_k;ubnH2_k;ubqH2Oloss_k;ubT_k;ubMbt;...
    ubT_bt_out;ubT_el_in;ubT_cw_out;ubU_el_k;ubq_lye_k;ubq_cw;ubqH2O]; 
 

%% preparing symbolic constraints
g = {};
% preparing numeric bounds
lbg = [];
ubg = [];

% declaring constraints

Iden = SX.zeros(par.N,1);
for nEl = 1:par.N
    Iden(nEl) = (0.1*xAlg(par.N+nEl))/par.EL(nEl).A; %current density in mA/cm2
end
IdenMin = 32;   %minimum current density, 32 mA/cm2
IdenMax = 198.5;%maximum current density, 198.5 mA/cm2

g = {g{:},eqnAlg, eqnDiff, Iden,eqnPnet};
lbg = [lbg;zeros(11*par.N,1);IdenMin*ones(par.N,1);0];
ubg = [ubg;zeros(11*par.N,1);IdenMax*ones(par.N,1);P0];


Objvol_H2 = SX.zeros(par.N,1);
for nEl = 1:par.N
    Objvol_H2(nEl) = (xAlg(4*par.N+nEl)*0.0224136*3600);%[Nm3/h]
end
J = -(Objvol_H2(1)+Objvol_H2(2)+Objvol_H2(3));

% J = 10;

%% formalize into an NLP problem
nlp = struct('x',vertcat(w{:}),'g',vertcat(g{:}),'f',J,'p',Pnet);


% assign solver - IPOPT in this case
solver = nlpsol('solver','ipopt',nlp);


% solve - using the defined initial guess and bounds
sol = solver('x0',w0,'lbx',lbw,'ubx',ubw,'lbg',lbg,'ubg',ubg,'p',P0);
res = full(sol.x);

%% Extracting results
Uk = [];
Ik = [];
Pk = [];
Feffk = [];
nH2k = [];
qH2Olossk = [];
Tk = [];
massBt=[];
T_bt_out=[];
T_el_in=[];        %temp of inlet lye stream coming into the electrolyzer
T_CW_out=[];
Vss = [];
q_lyek = [];
qf_cw = [];
Qwater = [];

for nEl = 1:par.N
    %optimal value of the algebriac state
    Uk = [Uk res(nEl)];                         %cell voltage of the electrolyzer
    Ik = [Ik res(par.N+nEl)];                   %current in the electrolyzer
    Pk = [Pk res(2*par.N+nEl)];                 %power of the individual electrolyzer
    Feffk = [Feffk res(3*par.N+nEl)];           %faraday efficiency of each electrolyzer
    nH2k = [nH2k res(4*par.N+nEl)];             %hydrogen produced form each individual electrolyzer
    qH2Olossk = [qH2Olossk res(5*par.N+nEl)];   %water loss during electrolysis in kth electrolyzer
    %optimal value of the differential state
    Tk = [Tk res(6*par.N+nEl)];               %temperature of the individual electrolyzer
    massBt = [massBt res(7*par.N+nEl)];
    T_bt_out = [T_bt_out res(8*par.N+nEl)];
    T_el_in = [T_el_in res(9*par.N+nEl)];         %temp of inlet lye stream coming into the electrolyzer
    T_CW_out = [T_CW_out res(10*par.N+nEl)];
    %optimal value of the inputs
    Vss = [Vss res(11*par.N+nEl)];            %electrolyzer voltage
    q_lyek = [q_lyek res(12*par.N+nEl)];      %lye flowrate
    qf_cw = [qf_cw res(13*par.N+nEl)];
    Qwater = [Qwater res(14*par.N+nEl)];
end

%% Calculation of initial state vector

%Nominal load H2 production and specific electricity consumption
V_H2_ini = nH2k*0.0224136*3600;%[Nm3/h]
for nEl = 1:par.N
    Ps_ini(nEl) = (Uk(nEl)*Ik(nEl)*par.EL(nEl).nc)/(1000*V_H2_ini(nEl));%[kWh/Nm3]
    Iden(nEl) = 0.1*Ik(nEl)/par.EL(nEl).A; 
end
Pnet = sum(Pk)

Iden
Tk
T_el_in

V_H2_ini
Ps_ini;

Eff_El = 3.55./Ps_ini



z0 = [Uk Ik Pk Feffk nH2k qH2Olossk];
x0 = [Tk massBt T_bt_out T_el_in T_CW_out];
u0 = [Vss q_lyek qf_cw Qwater];

end