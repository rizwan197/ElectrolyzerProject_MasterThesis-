function [z0, x0, u0, EXIT] = El_SteadyStateOptimization(N,X0,P0)

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
[lbz,lbx,lbu,ubz,ubx,ubu] = decision_varbound(par.N);
lbw = [lbw;lbz;lbx;lbu];
ubw = [ubw;ubz;ubx;ubu];


%% preparing symbolic constraints
g = {};
% preparing numeric bounds
lbg = [];
ubg = [];

% declaring additional constraints

Iden = SX.zeros(par.N,1);
for nEl = 1:par.N
    Iden(nEl) = (0.1*xAlg(par.N+nEl))/par.EL(nEl).A; %current density in mA/cm2
end
IdenMin = 32;   %minimum current density, 32 mA/cm2
IdenMax = 198.5;%maximum current density, 198.5 mA/cm2

deltaT1_k = xDiff(3*par.N+1:4*par.N) - par.Tw_in_k*ones(par.N,1);   %difference in temperatures of hot and cold streams at inlet
deltaT2_k = xDiff(2*par.N+1:3*par.N) - xDiff(4*par.N+1:5*par.N);    %difference in temperatures of hot and cold streams at outlet

deltaT_El1 = xDiff(1)-xDiff(3*par.N+1);
deltaT_El2 = xDiff(2)-xDiff(3*par.N+2);
deltaT_El3 = xDiff(3)-xDiff(4*par.N);

g = {g{:},eqnAlg, eqnDiff, Iden,eqnPnet,deltaT1_k,deltaT2_k,deltaT_El1,deltaT_El2,deltaT_El3};
lbg = [lbg;zeros(11*par.N,1);IdenMin*ones(par.N,1);0;2e-3*ones(par.N,1);2e-3*ones(par.N,1);0;0;0];
ubg = [ubg;zeros(11*par.N,1);IdenMax*ones(par.N,1);P0;inf*ones(par.N,1);inf*ones(par.N,1);30;30;30];

Objvol_H2 = SX.zeros(par.N,1);
for nEl = 1:par.N
    Objvol_H2(nEl) = (xAlg(4*par.N+nEl)*0.0224136*3600);%[Nm3/h]
end

J = -(Objvol_H2(1)+Objvol_H2(2)+Objvol_H2(3));
% J = -(Objvol_H2(1));

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

% Ps_ini; 
% Eff_El = 3.55./Ps_ini

z0 = [Uk Ik Pk Feffk nH2k qH2Olossk];
x0 = [Tk massBt T_bt_out T_el_in T_CW_out];
u0 = [Vss q_lyek qf_cw Qwater];

end