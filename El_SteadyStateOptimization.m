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


lbT_k = 25*ones(par.N,1);%lower bound on the electrolyzer temperature
ubT_k = 80*ones(par.N,1);

%constraints on the inputs
lbq_lye_k = 500*ones(par.N,1);%lower bound on the lye flowrate
ubq_lye_k = 8000*ones(par.N,1);

lbw = [lbw;lbu_k;lbi_k;lbP_k;lbFeff_k;lbnH2_k;lbT_k;lbq_lye_k];%bounds on all the variables
ubw = [ubw;ubu_k;ubi_k;ubP_k;ubFeff_k;ubnH2_k;ubT_k;ubq_lye_k]; 
 

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
lbg = [lbg;zeros(7*par.N+11,1);IdenMin*ones(par.N,1);0];
ubg = [ubg;zeros(7*par.N+11,1);IdenMax*ones(par.N,1);P0];


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

Tk = [];

q_lyek = [];

for nEl = 1:par.N
    %optimal value of the algebriac state
    Uk = [Uk res(nEl)];                         %cell voltage of the electrolyzer
    Ik = [Ik res(par.N+nEl)];                   %current in the electrolyzer
    Pk = [Pk res(2*par.N+nEl)];                 %power of the individual electrolyzer
    Feffk = [Feffk res(3*par.N+nEl)];           %faraday efficiency of each electrolyzer
    nH2k = [nH2k res(4*par.N+nEl)];             %hydrogen produced form each individual electrolyzer
    
    %optimal value of the differential state
    Tk = [Tk res(5*par.N+nEl)];               %temperature of the individual electrolyzer
    
    %optimal value of the inputs
    q_lyek = [q_lyek res(6*par.N+nEl)];      %lye flowrate
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



z0 = [Uk Ik Pk Feffk nH2k];
x0 = [Tk];
u0 = [q_lyek];

end