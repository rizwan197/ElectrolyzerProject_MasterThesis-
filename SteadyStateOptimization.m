%This file solves the steady state solution for electrolyzer plant
%All the dynamic simulation states are solved at steady state

%nEl = sequence of the electrolyzer

clc
clear

N=3;
par = parElectrolyzer(N);


%% Load CasADi
import casadi.*

%% Define symbolic variables
x = SX.sym('x',7*par.N+10);              %symbolic variables for cell voltage, current and electrolyzer voltage (V)
eqnAlg = SX.zeros(6*par.N+5,1);
eqnDiff = SX.zeros(par.N+5,1);

%variables for algebriac eqns.
u_k=[];
i_k=[];
P_k=[];
Feff_k=[];
nH2_k=[];
qH2Oloss_k=[];
nH2El_net=x(6*par.N+1);
nH2out_net=x(6*par.N+2);
nO2El_net=x(6*par.N+3);
nO2out_net=x(6*par.N+4);
T_El_out=x(6*par.N+5);                  %temp after mixing of the exiting liquid streams from all electrolyzers

%variables for differential equations 
T_k=[];
Psto_H2=x(7*par.N+6);
Psto_O2=x(7*par.N+7);
Mass_Bt=x(7*par.N+8);
T_El_in=x(7*par.N+9);                   %temp of inlet lye stream coming into the electrolyzer
T_cw_out=x(7*par.N+10);


for nEl = 1:par.N
    u_k = [u_k x(nEl)];                         %cell voltage of the electrolyzer
    i_k = [i_k x(par.N+nEl)];                   %current in the electrolyzer
    P_k = [P_k x(2*par.N+nEl)];                 %power of the individual electrolyzer
    Feff_k = [Feff_k x(3*par.N+nEl)];           %faraday efficiency of each electrolyzer
    nH2_k = [nH2_k x(4*par.N+nEl)];             %hydrogen produced form each individual electrolyzer
    qH2Oloss_k = [qH2Oloss_k x(5*par.N+nEl)];   %water loss during electrolysis in kth electrolyzer
    T_k = [T_k x(6*par.N+5+nEl)];               %temperature of the individual electrolyzer
end
     

%% Define inputs for the simulation (MVs for the dynamic simulation)

inp = SX.sym('inp',2*par.N+4);

U_El_k=[];
q_lye_k=[];
for nEl = 1:par.N
    U_El_k = [U_El_k inp(nEl)];
    q_lye_k = [q_lye_k inp(par.N+nEl)];
end

q_cw = inp(2*par.N+1);
zH2 = inp(2*par.N+2);
zO2 = inp(2*par.N+3);
q_H2O = inp(2*par.N+4);
 
%% Initial guess for steady state calculations
%algebriac state variables
u_k0 = 1.8*ones(1,par.N);       %initial guess for cell voltage
P_k0 = 2135000*ones(1,par.N);
i_k0 = P_k0./(u_k0.*par.EL(1).nc);     %initial guess for current
Feff_k0 = 0.97*ones(1,par.N);
nH2_k0 = 6*ones(1,par.N);       %[mol/s]
qH2Oloss_k0 = nH2_k0*par.Const.Mwt.*ones(1,par.N);%[g/s]
nH2El_net0 = sum(nH2_k0);       %[mol/s]
nH2out_net0 = sum(nH2_k0);
nO2El_net0 = 0.5*nH2El_net0;
nO2out_net0 = 0.5*nH2out_net0;
T_El_out0 = 80;                 %initial guess for the temperature of lye entering in the heat exchanger

%differential state variables
T_k0 = 75*ones(1,par.N);
Psto_H20 = 20;      %initial H2 storage pressure (calculated from steady state solution) [bar]
Psto_O20 = 20;      %initial O2 storage pressure (calculated from steady state solution) [bar]
Mass_Bt0 = 1000000; %mass of the water in the buffer tank,[g]
T_El_in0 = 65;      %initial guess for the temperature of inlet lye into the electrolyzer, [deg C]
T_cw_out0 = 20;     %initial guess for the exit temperature of the cooling water leaving heat exchanger, [deg C]

z0 = [u_k0 i_k0 P_k0 Feff_k0 nH2_k0 qH2Oloss_k0 nH2El_net0 nH2out_net0 nO2El_net0 nO2out_net0 T_El_out0];
x0 = [T_k0 Psto_H20 Psto_O20 Mass_Bt0 T_El_in0 T_cw_out0];

%initial guess for input variables
U_El_k_0 = 414.0301*ones(1,par.N);    %voltage across electrolyzers, [Volts]
q_lye_k_0 = 6648*ones(1,par.N);
q_cw_0 = 2.0698e4;                    %cooling water flow rate, [g/s]
zH2_0 = 0.4;
zO2_0 = 0.4;
q_H2O_0 = 324.2657;                    %total water lost during electrolysis, [grams/sec]

Par = [U_El_k_0 q_lye_k_0 q_cw_0 zH2_0 zO2_0 q_H2O_0]; 

%initial guess vector for the IPOPT
X0 = [z0 x0 Par];

%% Model equations
for nEl = 1:par.N
    
    %relation for Urev with T from LeRoy eqn. 58
    par.EL(nEl).Urev = 1.5184 - 1.5421e-3*(273+T_k(nEl)) + 9.523e-5*(273+T_k(nEl))*log((273+T_k(nEl))) + ...
        9.84e-8*(273+T_k(nEl))^2;
    
    %model equations
    eqnAlg(nEl) = u_k(nEl)*i_k(nEl)*par.EL(nEl).nc - P_k(nEl);                        %power = nc*UI
    eqnAlg(par.N+nEl) = u_k(nEl) - (par.U(nEl).r1+par.U(nEl).r2*T_k(nEl))*i_k(nEl)/par.EL(nEl).A - par.U(nEl).s*log10(((par.U(nEl).t1+par.U(nEl).t2/T_k(nEl)+...
        par.U(nEl).t3/(T_k(nEl)^2))*i_k(nEl)/par.EL(nEl).A)+1) - par.EL(nEl).Urev;  %U-I relationship
    eqnAlg(2*par.N+nEl) = u_k(nEl)*par.EL(nEl).nc - U_El_k(nEl);                    %U(j).nc(j)=V
    eqnAlg(3*par.N+nEl) = Feff_k(nEl) - ((.1*i_k(nEl)/par.EL(nEl).A)^2)/(par.U(nEl).f1+((.1*i_k(nEl)/par.EL(nEl).A)^2))*par.U(nEl).f2;     %faraday efficiency
    eqnAlg(4*par.N+nEl) = nH2_k(nEl) - Feff_k(nEl)*par.EL(nEl).nc*i_k(nEl)/(par.Const.ze*par.Const.FC);%nH2, H2 production rate from individual electrolyzer
    eqnAlg(5*par.N+nEl) = qH2Oloss_k(nEl) - nH2_k(nEl)*par.Const.Mwt;  %flowrate of water lost, [g/s]  
end

sum_H2net = SX.zeros(1,1);
netqlye = SX.zeros(1,1);
netqloss = SX.zeros(1,1);
qlyeT = SX.zeros(1,1);
qlossT = SX.zeros(1,1);

for nEl = 1:par.N
    sum_H2net = sum_H2net + nH2_k(nEl);       %sum of hydrogen from all individual electrolyzers
    netqlye = netqlye + q_lye_k(nEl);             %sum the lye flowing into all individual electrolyzers
    netqloss = netqloss + qH2Oloss_k(nEl);           %sum of water lost from all individual electrolyzers
    qlyeT = qlyeT + q_lye_k(nEl)*T_k(nEl);              %calculate term qlye(k).T(k)
    qlossT = qlossT + qH2Oloss_k(nEl)*T_k(nEl);            %calculate term qH2Oloss(k).T(k)
end


eqnAlg(6*par.N+1) = nH2El_net - sum_H2net;                                              %algebraic eqn for net hydrogen flowrate from all the electrolyzers
eqnAlg(6*par.N+2) = nH2out_net - par.kvalveH2*zH2*sqrt(Psto_H2-par.Storage.PoutH2);     %algebraic eqn for net hydrogen flowrate from the storage tank
eqnAlg(6*par.N+3) = nO2El_net - nH2El_net/2;                                            %algebraic eqn for net oxygen flowrate from all the electrolyzers
eqnAlg(6*par.N+4) = nO2out_net - par.kvalveO2*zO2*sqrt(Psto_O2-par.Storage.PoutO2);     %algebraic eqn for net oxygen flowrate from the storage tank
eqnAlg(6*par.N+5) = T_El_out - (((qlyeT*par.Const.CpLye) - (qlossT*par.Const.Cp) + netqloss*(par.Const.Cp-par.Const.CpLye)*par.Const.Tref)/...
    ((netqlye-netqloss)*par.Const.CpLye));                                          %calculation of Tout

for nEl = 1:par.N
    eqnDiff(nEl) = (q_lye_k(nEl)*par.Const.CpLye*(T_El_in-T_k(nEl)) + par.EL(nEl).nc*(u_k(nEl)-par.EL(nEl).Utn)*i_k(nEl) - ...
        par.TherMo(nEl).A_surf*i_k(nEl)* par.EL(nEl).nc/1000*(par.TherMo(nEl).hc*(T_k(nEl)-par.EL(nEl).Ta) + ...
        par.sigma*par.em*((T_k(nEl)+273.15)^4-(par.EL(nEl).Ta+273.15)^4)))/(par.TherMo(nEl).CtS*i_k(nEl));%differential eqn for the electrolyzer temperature
end

eqnDiff(par.N+1) = (par.Storage.TstoH2*par.Storage.Rg/par.Storage.VstoH2)*(nH2El_net-nH2out_net);  %differential eqn for hydrogen storage pressure
eqnDiff(par.N+2) = (par.Storage.TstoO2*par.Storage.Rg/par.Storage.VstoO2)*(nO2El_net-nO2out_net);  %differential eqn for oxygen storage pressure
eqnDiff(par.N+3) = (netqlye-netqloss) + q_H2O - netqlye;                                     %differential eqn for mass in the buffer tank, [grams i.e. pho*V]

%dynamic thermal balance for the heat exchanger, T_in and Tw_out are differential variables

deltaT1 = T_El_in - par.Tw_in;%difference in temperatures of hot and cold streams at inlet
deltaT2 = T_El_out - T_cw_out;%difference in temperatures of hot and cold streams at outlet
deltaT_LMTD = (deltaT1-deltaT2)/log(deltaT1/deltaT2);

eqnDiff(par.N+4) = ((netqlye*(T_El_out-T_El_in))/(1000*par.Const.rho*par.Const.Vh)) - (par.Hex.UA*deltaT_LMTD/(1000*par.Const.rho*par.Const.CpLye*par.Const.Vh));%differential eqn for the hot stream exit temp from heat exchanger
eqnDiff(par.N+5) = ((q_cw*(par.Tw_in-T_cw_out))/(1000*par.Const.rho*par.Const.Vc)) + (par.Hex.UA*deltaT_LMTD/(1000*par.Const.rho*par.Const.Cp*par.Const.Vc));%differential eqn for the cold stream exit temperature from heat exchanger

%% Solving steady state using rootfinder (Can't be used to solve steady state optimimization problem, since we have inputs as well)

% g = Function('g',{x},{eqnAlg,eqnDiff});
% G = rootfinder('G','newton',g);
% res = full(G(X0));


%% Solving steady state using IPOPT (constant value as objective)
%Nominal load H2 production and specific electricity consumption
% V_H2_ini = SX.zeros(1,par.N);
% Ps_ini = SX.zeros(1,par.N);
% 
% for nEl = 1:par.N
%     Ps_ini(nEl) = (u_k(nEl)*i_k(nEl)*par.EL(nEl).nc)/(1000*V_H2_ini(nEl));%[kWh/Nm3]
%     V_H2_ini(nEl) = nH2_k(nEl)*0.0224136*3600;%[Nm3/h]
% end
% L=((sum(V_H2_ini)-485)^2);

L=1;

% preparing symbolic variables
w = {};
% preparing numeric variables and bounds
w0 = [];
lbw = [];
ubw = [];

% declaring them symbolic
w = {w{:},x,inp};
% declaring them numerically
w0 = [w0;X0']; %initial guess
%defining constraints on the decision variables (states and inputs i.e. MVs)
%constraints on states
lbu_k = zeros(par.N,1);%lower bound on cell voltage
ubu_k = 1.9*ones(par.N,1);
lbi_k = zeros(par.N,1);%lower bound on current
ubi_k = inf*ones(par.N,1);
lbP_k = zeros(par.N,1);%lower bound on power
ubP_k = inf*ones(par.N,1);
lbFeff_k = zeros(par.N,1);%lower bound on faraday efficiency
ubFeff_k = inf*ones(par.N,1);
lbnH2_k = zeros(par.N,1);%lower bound on hydrogen production
ubnH2_k = inf*ones(par.N,1);
lbqH2Oloss_k = zeros(par.N,1);%lower bound on water loss
ubqH2Oloss_k = inf*ones(par.N,1);
lbnH2el_net = 0;%lower bound on net hydrogen production from the electrolyzer
ubnH2el_net = inf;
lbnH2out_net = 0;%lower bound on hydrogen from the outlet 
ubnH2out_net = inf;
lbnO2el_net = 0;%lower bound on net oxygen production from the electrolyzer
ubnO2el_net = inf;
lbnO2out_net = 0;%lower bound on oxygen from the outlet
ubnO2out_net = inf;
lbT_el_out = 0;%lower bound on the mixed lye stream temperature at electrolyzer outlet
ubT_el_out = inf;
lbT_k = 65*ones(par.N,1);%lower bound on the electrolyzer temperature
ubT_k = 90*ones(par.N,1);
lbPstoH2 = 0;%lower bound on the hydrogen storage pressure
ubPstoH2 = inf;
lbPstoO2 = 0;%lower bound on the oxygen storage pressure
ubPstoO2 = inf;
lbMbt = 0;%lower bound on the mass in the buffer tank 
ubMbt = inf;
lbT_el_in = 0;%lower bound on the mixed lye stream temperature at electrolyzer inlet
ubT_el_in = inf;
lbT_cw_out = 0;%lower bound on the coolant outlet temperature
ubT_cw_out = inf;

% z0 = [u_k0 i_k0 P_k0 Feff_k0 nH2_k0 qH2Oloss_k0 nH2El_net0 nH2out_net0 nO2El_net0 nO2out_net0 T_El_out0];
% x0 = [T_k0 Psto_H20 Psto_O20 Mass_Bt0 T_El_in0 T_cw_out0];

%constraints on the inputs
lbU_el_k = zeros(par.N,1);%lower bound on the electrolyzer voltage
ubU_el_k = inf*ones(par.N,1);
lbq_lye_k = zeros(par.N,1);%lower bound on the lye flowrate
ubq_lye_k = inf*ones(par.N,1);
lbq_cw = 0;%lower bound on the coolant flow rate
ubq_cw = inf;
lbzH2 = 0;%lower bound on hydrogen outlet valve opening
ubzH2 = inf;
lbzO2 = 0;%lower bound on oxygen outlet valve opening
ubzO2 = inf;
lbqH2O = 0;%lower bound on total water lost during electrolysis
ubqH2O = inf;

lbw = [lbw;lbu_k;lbi_k;lbP_k;lbFeff_k;lbnH2_k;lbqH2Oloss_k;...
    lbnH2el_net;lbnH2out_net;lbnO2el_net;lbnO2out_net;lbT_el_out;lbT_k;lbPstoH2;lbPstoO2;lbMbt;lbT_el_in;lbT_cw_out;...
    lbU_el_k;lbq_lye_k;lbq_cw;lbzH2;lbzO2;lbqH2O];%bounds on all the variables
ubw = [ubw;ubu_k;ubi_k;ubP_k;ubFeff_k;ubnH2_k;ubqH2Oloss_k;...
    ubnH2el_net;ubnH2out_net;ubnO2el_net;ubnO2out_net;ubT_el_out;ubT_k;ubPstoH2;ubPstoO2;ubMbt;ubT_el_in;ubT_cw_out;...
    ubU_el_k;ubq_lye_k;ubq_cw;ubzH2;ubzO2;ubqH2O]; 
 

% preparing symbolic constraints
g = {};
% preparing numeric bounds
lbg = [];
ubg = [];

% declaring constraints
g = {g{:},eqnAlg, eqnDiff};
lbg = [lbg;zeros(7*par.N+10,1)];
ubg = [ubg;zeros(7*par.N+10,1)];

% optimization objective function
% By default, casadi always minimizes the problem. 
% Since we want to minimize it, we have to write:
J = ([x;inp]-w0)'*([x;inp]-w0); %since we want to find optimal near the initial guess
% J = L;

% formalize it into an NLP problem
% the main difference is the fourth argument 'p'!
nlp = struct('x',vertcat(w{:}),'g',vertcat(g{:}),'f',J);

% assign solver - IPOPT in this case
solver = nlpsol('solver','ipopt',nlp);

% solve - using the previous defined initial guess and bounds
sol = solver('x0',w0,'lbx',lbw,'ubx',ubw,'lbg',lbg,'ubg',ubg);
res = full(sol.x);

Uk = [];
Ik = [];
Pk = [];
Feffk = [];
nH2k = [];
qH2Olossk = [];
Tk = [];

for nEl = 1:par.N
    Uk = [Uk res(nEl)];                         %cell voltage of the electrolyzer
    Ik = [Ik res(par.N+nEl)];                   %current in the electrolyzer
    Pk = [Pk res(2*par.N+nEl)];                 %power of the individual electrolyzer
    Feffk = [Feffk res(3*par.N+nEl)];           %faraday efficiency of each electrolyzer
    nH2k = [nH2k res(4*par.N+nEl)];             %hydrogen produced form each individual electrolyzer
    qH2Olossk = [qH2Olossk res(5*par.N+nEl)];   %water loss during electrolysis in kth electrolyzer
    Tk = [Tk res(6*par.N+5+nEl)];               %temperature of the individual electrolyzer
end

nH2El_tot=res(6*par.N+1);
nH2out_tot=res(6*par.N+2);
nO2El_tot=res(6*par.N+3);
nO2out_tot=res(6*par.N+4);
T_el_out=res(6*par.N+5);        %temp after mixing of the exiting liquid streams from all electrolyzers
PstoH2=res(7*par.N+6);
PstoO2=res(7*par.N+7);
massBt=res(7*par.N+8);
T_el_in=res(7*par.N+9);         %temp of inlet lye stream coming into the electrolyzer
T_CW_out=res(7*par.N+10);


%% Calculation of initial state vector

%Nominal load H2 production and specific electricity consumption
V_H2_ini = nH2k*0.0224136*3600;%[Nm3/h]
for nEl = 1:par.N
    Ps_ini(nEl) = (Uk(nEl)*Ik(nEl)*par.EL(nEl).nc)/(1000*V_H2_ini(nEl));%[kWh/Nm3]
end

z_ss = [Uk Ik Pk Feffk nH2k qH2Olossk nH2El_tot nH2out_tot nO2El_tot nO2out_tot T_el_out];
x_ss = [Tk PstoH2 PstoO2 massBt T_el_in T_CW_out];

for nEl = 1:par.N
    V_H2_ini(nEl) = nH2k(nEl)*0.0224136*3600;%[Nm3/h]
    Ps_ini(nEl) = (Uk(nEl)*Ik(nEl)*par.EL(nEl).nc)/(1000*V_H2_ini(nEl));%[kWh/Nm3]
end
V_H2_ini
Ps_ini;