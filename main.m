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
T_bt_in0 = 73*ones(1,par.N);            %Initial guess for the temperature of lye mixture at the inlet of the buffer tank,[degC]

%differential state variables('x')
T_k0 = 75*ones(1,par.N);
Mass_Bt0 = 2000000*ones(1,par.N);  %mass of the liquid in the buffer tank,[g] 2000kg/Electrolyzer
T_bt_out0 = 70*ones(1,par.N);       %Initial guess for the temperature of lye mixture at the exit of the buffer tank,[degC]
T_El_in0 = 65*ones(1,par.N);        %initial guess for the temperature of inlet lye into the electrolyzer, [deg C]
T_cw_out0 = 20*ones(1,par.N);       %initial guess for the exit temperature of the cooling water leaving heat exchanger,[deg C]

z_guess = [u_k0 i_k0 P_k0 Feff_k0 nH2_k0 qH2Oloss_k0 T_bt_in0];
x_guess = [T_k0 Mass_Bt0 T_bt_out0 T_El_in0 T_cw_out0];

%initial guess for input variables('u')
U_El_k_0 = 414.0301*ones(1,par.N);      %voltage across electrolyzers, [Volts]
q_lye_k_0 = 6648*ones(1,par.N);         %lye flowrate, [g/s]
q_cw_0 = 2.0698e4/N*ones(1,par.N);                      %cooling water flow rate, [g/s]
q_H2O_0 = 324.2657/N*ones(1,par.N);                     %total water lost during electrolysis, [grams/sec]

u_guess = [U_El_k_0 q_lye_k_0 q_cw_0 q_H2O_0]; 
counter = 1;
flag = {};

for Pnet = 9e6:-0.1e6:.6e6
    
%initial guess vector for the IPOPT
X_guess = [z_guess x_guess u_guess];

%% Solve the steady state problem
[z0, x0, u0, EXIT] = El_SteadyStateOptimization(N,X_guess,Pnet);
z_guess = z0;
x_guess = x0;
u_guess = u0;

T_El_in_k = x0(3*par.N+1:4*par.N);%temperature of lye entering the electrolyzer
T_bt_out_k = x0(2*par.N+1:3*par.N);
T_cw_out_k = x0(4*par.N+1:5*par.N);

%Initial value of the MVs 
Vss = u0(1:par.N);
q_lyek = u0(par.N+1:2*par.N);
qlye_kgs = q_lyek/1000;
qf_cw = u0(2*par.N+1:3*par.N);
qcw_kgs = qf_cw/1000;
Qwater = u0(3*par.N+1:4*par.N);

Pcons = sum(z0(2*par.N+1:3*par.N));
Iden = (0.1*z0(par.N+1:2*par.N))/par.EL(1).A;
Tk = x0(1:par.N);
V_H2_ini = z0(4*par.N+1:5*par.N)*0.0224136*3600;


if strcmp(EXIT, 'Solve_Succeeded')
    ac_DC_S2(counter,:) = [Iden/198.5, 32./Iden, Tk/80, 25./Tk, (Tk-T_El_in_k)/30, 2e-3./(T_El_in_k - par.Tw_in_k*ones(1,par.N)), 2e-3./(T_bt_out_k-T_cw_out_k),...
        qlye_kgs/10, 0.5./qlye_kgs, qcw_kgs/20 1e-5./qcw_kgs];
    row_DC_S2_New(counter,:) = [Pnet/1e6,Pcons/1e6,qlye_kgs,qcw_kgs,Iden,Tk,T_El_in_k,T_cw_out_k,T_bt_out_k,V_H2_ini,sum(V_H2_ini)];
    
else
    ac_DC_S2(counter,:) = NaN*ones(1,11*par.N);
    row_DC_S2_New(counter,:) = 0*ones(1,8*par.N+3);
end

flag = {flag{:},EXIT}';
counter = counter+1;
end

% save('Data_DCEl_S2_NewHex')

