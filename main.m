clc
clear
close all

%% Load CasADi
addpath('/Users/mdrizwan/Documents/MATLAB/casadi-osx-matlabR2015a-v3.5.3')
import casadi.*

%% Loading parameters
N = 3;                               %no. of electrolyzers
par = parElectrolyzer(N);

%% Inputs for the simulation
num_hr = .25;                           %no. of hours
t0 = 1;                                 %start, [s)]
ts = 1;                                 %time step, [s]
tf = num_hr*60*60;                      %final, [s]
tsamp = t0:ts:tf;
len = length(tsamp);                    %number of simulation time steps
tstep = 200;

%% Initial guess for steady state solution using IPOPT

%disturbance is total power
Pnet = 9e6;%1900e3*par.N; %6.3MW total input power

%algebriac state variables('z')
u_k0 = 1.8*ones(1,par.N);               %initial guess for cell voltage
P_k0 = Pnet/par.N*ones(1,par.N);        %intial guess is power divided equally among electrolyzers 
i_k0 = P_k0./(u_k0.*par.EL(1).nc);      %initial guess for current
Feff_k0 = 0.97*ones(1,par.N);
nH2_k0 = 6*ones(1,par.N);               %[mol/s]
qH2Oloss_k0 = nH2_k0*par.Const.Mwt.*ones(1,par.N);%[g/s]
nH2El_net0 = sum(nH2_k0);               %[mol/s]
nH2out_net0 = sum(nH2_k0);
nO2El_net0 = 0.5*nH2El_net0;
nO2out_net0 = 0.5*nH2out_net0;
T_El_out0 = 80;                         %initial guess for the temperature of lye entering in the heat exchanger

%differential state variables('x')
T_k0 = 75*ones(1,par.N);
Psto_H20 = 25;        %initial H2 storage pressure (calculated from steady state solution) [bar]
Psto_O20 = 25;        %initial O2 storage pressure (calculated from steady state solution) [bar]
Mass_Bt0 = 6000000;  %mass of the liquid in the buffer tank,[g] 6000kg
T_bt_out0 = 70;       %Initial guess for the temperature of lye mixture at the exit of the buffer tank,[degC]
T_El_in0 = 65;        %initial guess for the temperature of inlet lye into the electrolyzer, [deg C]
T_cw_out0 = 20;       %initial guess for the exit temperature of the cooling water leaving heat exchanger,[deg C]

z_guess = [u_k0 i_k0 P_k0 Feff_k0 nH2_k0 qH2Oloss_k0 nH2El_net0 nH2out_net0 nO2El_net0 nO2out_net0 T_El_out0];
x_guess = [T_k0 Psto_H20 Psto_O20 Mass_Bt0 T_bt_out0 T_El_in0 T_cw_out0];

%initial guess for input variables('u')
U_El_k_0 = 414.0301*ones(1,par.N);      %voltage across electrolyzers, [Volts]
q_lye_k_0 = 6648*ones(1,par.N);         %lye flowrate, [g/s]
q_cw_0 = 2.0698e4;                      %cooling water flow rate, [g/s]
zH2_0 = 0.4;
zO2_0 = 0.4;
q_H2O_0 = 324.2657;                     %total water lost during electrolysis, [grams/sec]

u_guess = [U_El_k_0 q_lye_k_0 q_cw_0 zH2_0 zO2_0 q_H2O_0]; 

counter = 1;
flag = {};

for Pnet = 9e6:-0.1e6:.6e6
    
X_guess = [z_guess x_guess u_guess];

%% Solve the steady state problem
[z0, x0, u0, EXIT] = El_SteadyStateOptimization(N,X_guess,Pnet);

z_guess = z0;
x_guess = x0;
u_guess = u0;

T_El_in_set = x0(par.N+5);%setpoint for the temperature of lye entering the electrolyzer 
T_cw_out = x0(par.N+6);
T_bt_out = x0(par.N+4);

%Initial value of the MVs 
Vss = u0(1:par.N);
q_lyek = u0(par.N+1:2*par.N);
qlye_kgs = q_lyek/1000; 
qf_cw = u0(2*par.N+1);
qcw_kgs = qf_cw/1000;
zH2 = u0(2*par.N+2);
zO2 = u0(2*par.N+3);
Qwater = u0(2*par.N+4);
nH2 = sum(z0(4*par.N+1:5*par.N));

Pcons = sum(z0(2*par.N+1:3*par.N));
Iden = 0.1*z0(par.N+1:2*par.N)./par.EL(1).A;
Tk = x0(1:par.N);
V_H2_ini = z0(4*par.N+1:5*par.N)*0.0224136*3600;

for nEl = 1:par.N
    Qlyeloss(nEl) = (q_lyek(nEl)*par.Const.CpLye*(T_El_in_set-Tk(nEl)));
    Qgenk(nEl) = par.EL(nEl).nc*(z0(nEl)-par.EL(nEl).Utn)*z0(par.N+nEl);
    Qlossk(nEl) = par.TherMo(nEl).A_El*(par.TherMo(nEl).hc*(Tk(nEl)-par.EL(nEl).Ta) + par.sigma*par.em*((Tk(nEl)+273.15)^4-(par.EL(nEl).Ta+273.15)^4));
    Qnet(nEl) = Qlyeloss(nEl)+Qgenk(nEl)-Qlossk(nEl);
end

row_C_S2(counter,:) = [Pnet/1e6,Pcons/1e6,qlye_kgs,qcw_kgs,Iden,Tk,T_El_in_set,T_cw_out,T_bt_out,V_H2_ini, sum(V_H2_ini) Qnet];

if strcmp(EXIT,'Solve_Succeeded')
    ac_C_S2(counter,:) = [Iden/198.5, 32./Iden, Tk/80, 25./Tk, (max(Tk)-T_El_in_set)/30, Pcons/Pnet,...
        qlye_kgs/10, 0.5./qlye_kgs, qcw_kgs/80 1e-5/qcw_kgs];
    plot_C_S2_DegHex(counter,:) = [Pnet/1e6,Pcons/1e6,qlye_kgs,qcw_kgs,Vss,Iden,Tk,T_El_in_set,T_cw_out,T_bt_out,V_H2_ini, sum(V_H2_ini)];
else
    ac_C_S2(counter,:) = NaN*ones(1,6*par.N+4);
    plot_C_S2_DegHex(counter,:) = 0*ones(1,5*par.N+7);
end

flag = {flag{:},EXIT}';
counter = counter+1;

end

% save('DataThesis_ActiveRegionPlot')