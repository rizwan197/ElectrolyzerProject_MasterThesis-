function [z0, x0, u0] = El_RTO(N,X0,P0)

import casadi.*
par = parElectrolyzer(N);

%% Build the plant model and solve steady state optimization problem
[xDiff, xAlg, input, eqnAlg, eqnDiff] = modelnew(par.N);
x = [xAlg;xDiff]; 

%% Defining the disturbance
Pnet = MX.sym('Pnet');
Ptot = MX.zeros(1,1);

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
[lbx,lbz,lbu,ubx,ubz,ubu] = decision_var_bounds(par.N);

lbw = [lbw;lbz;lbx;lbu];
ubw = [ubw;ubz;ubx;ubu]; 

%% preparing symbolic constraints
g = {};
% preparing numeric bounds
lbg = [];
ubg = [];

% declaring additional constraints
uElconst = [];
for nEl=1:par.N-1
    uElconst = [uElconst;xAlg(nEl)-xAlg(nEl+1)];
end

deltaT1 = xDiff(par.N+5)-par.Tw_in;
deltaT2 = xDiff(par.N+4)-xDiff(par.N+6);

Iden = MX.zeros(par.N,1);
for nEl = 1:par.N
    Iden(nEl) = (0.1*xAlg(par.N+nEl))/par.EL(nEl).A; %current density in mA/cm2
end
IdenMin = 32;   %minimum current density, 32 mA/cm2
IdenMax = 198.5;%maximum current density, 198.5 mA/cm2

g = {g{:},eqnAlg, eqnDiff,uElconst, Iden,eqnPnet};
lbg = [lbg;zeros(7*par.N+11,1);zeros(par.N-1,1);IdenMin*ones(par.N,1);0];
ubg = [ubg;zeros(7*par.N+11,1);zeros(par.N-1,1);IdenMax*ones(par.N,1);P0];


Objvol_H2 = MX.zeros(par.N,1);
for nEl = 1:par.N
    Objvol_H2(nEl) = (xAlg(4*par.N+nEl)*0.0224136*3600);%[Nm3/h]
end
J = -(Objvol_H2(1)+Objvol_H2(2)+Objvol_H2(3));

% J = 10;

%% formalize into an NLP problem
nlpRTO = struct('x',vertcat(w{:}),'g',vertcat(g{:}),'f',J,'p',Pnet);

% assign solver - IPOPT in this case
solver = nlpsol('solver','ipopt',nlpRTO);

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