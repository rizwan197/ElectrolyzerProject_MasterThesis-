clc
clear
close all

%% Load CasADi
addpath('/Users/mdrizwan/Documents/MATLAB/casadi-osx-matlabR2015a-v3.5.1')
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
T_bt_in0 = 77*ones(1,par.N);       %Initial guess for the temperature of lye mixture at the inlet of the buffer tank,[degC]
T_El_in0 = 65*ones(1,par.N);        %initial guess for the temperature of inlet lye into the electrolyzer, [deg C]
T_cw_out0 = 12*ones(1,par.N);       %initial guess for the exit temperature of the cooling water leaving heat exchanger,[deg C]

%differential state variables('x')
T_k0 = 78*ones(1,par.N);
Mass_Bt0 = 2000000*ones(1,par.N);  %mass of the liquid in the buffer tank,[g] 2000kg/Electrolyzer
T_bt_out0 = 79*ones(1,par.N);       %Initial guess for the temperature of lye mixture at the exit of the buffer tank,[degC]


z_guess = [u_k0 i_k0 P_k0 Feff_k0 nH2_k0 qH2Oloss_k0 T_bt_in0 T_El_in0 T_cw_out0];
x_guess = [T_k0 Mass_Bt0 T_bt_out0];

%initial guess for input variables('u')
U_El_k_0 = 414.0301*ones(1,par.N);      %voltage across electrolyzers, [Volts]
q_lye_k_0 = 6648*ones(1,par.N);         %lye flowrate, [g/s]
% q_cw_0 = 2.0698e4/N*ones(1,par.N);                    %cooling water flow rate, [g/s]
Qcool_k_0 = ((8e4/par.N)*par.Const.Cp*(T_cw_out0 - 10));  %cooling duty,[J/s]
q_H2O_0 = 324.2657/N*ones(1,par.N);                     %total water lost during electrolysis, [grams/sec]

u_guess = [U_El_k_0 q_lye_k_0 Qcool_k_0 q_H2O_0]; 
counter = 1;
flag = {};

for Pnet = 9e6:-0.1e6:.6e6
    
%initial guess vector for the IPOPT
X_guess = [z_guess x_guess u_guess];

%% Solve the steady state problem
[z0, x0, u0, EXIT] = El_SteadyStateOptimizationEdit(N,X_guess,Pnet);
z_guess = z0;
x_guess = x0;
u_guess = u0;

T_bt_out_k = x0(2*par.N+1:3*par.N);
T_El_in_k = z0(7*par.N+1:8*par.N);%temperature of lye entering the electrolyzer
T_cw_out_k = z0(8*par.N+1:9*par.N);

%Initial value of the MVs 
Vss = u0(1:par.N);
q_lyek = u0(par.N+1:2*par.N);
qlye_kgs = q_lyek/1000;
Q_coolk = u0(2*par.N+1:3*par.N);%[J/s]
Q_coolk_kJs = Q_coolk/1e3;%[kJ/s]
Qwater = u0(3*par.N+1:4*par.N);

Pcons = sum(z0(2*par.N+1:3*par.N));
Iden = (0.1*z0(par.N+1:2*par.N))/par.EL(1).A;
Tk = x0(1:par.N);
V_H2_ini = z0(4*par.N+1:5*par.N)*0.0224136*3600;


if strcmp(EXIT, 'Solve_Succeeded')
%     ac_DC_S2(counter,:) = [Iden/198.5, 32./Iden, Tk/80, 25./Tk, (Tk-T_El_in_k)/30, 2e-3./(T_El_in_k - par.Tw_in_k*ones(1,par.N)), 2e-3./(T_bt_out_k-T_cw_out_k),...
%         qlye_kgs/10, 0.5./qlye_kgs, qcw_kgs/20 1e-5./qcw_kgs];
    row_DC_S2_New(counter,:) = [Pnet/1e6,Pcons/1e6,qlye_kgs,Q_coolk_kJs,Iden,Tk,T_El_in_k,T_cw_out_k,T_bt_out_k,V_H2_ini,sum(V_H2_ini)];
    
else
%     ac_DC_S2(counter,:) = NaN*ones(1,11*par.N);
    row_DC_S2_New(counter,:) = 0*ones(1,8*par.N+3);
end

flag = {flag{:},EXIT}';
counter = counter+1;
end
(8e4/N)*par.Const.Cp*70
% save('Data_DCEl_S2_NewHex')

%% Build the plant model
[xDiff, xAlg, input, eqnAlg, eqnDiff, F] = model(par.N);

%% Manipulated variables
%these are the degree of freedoms that we will utilise to control the system

V_El = zeros(len,N);              %voltage across the electrolyzer cell, [volts], len is the length of time vector
for j = 1:N
    V_El(1:end,j) = Vss(j)*1;     %incremental step change in common voltage across all electrolysers
    %     V_El(tstep:end,j)=Vss(j)*1;
end

qlye = zeros(len,N);                   %lye flowrate, [g/s]
for j = 1:N
    qlye(1:end,j) = q_lyek(j)*1;       %assumed same lye flowarate to all the electrolyzers
end
% qlye(tstep:end,2) = q_lyek(2)*1.2;

q_cw = zeros(len,N);
for j=1:N
q_cw(1:end,j) = qf_cw(j)*1;                     %cooling water flow rate as a manipulated variable, [g/s]
end

qH2O = zeros(len,N);
for j=1:N
qH2O(1:end,j) = Qwater(j)*1;  %flow rate of water added to buffer tank as a manipulated variable, [g/s]
end

%% Initialize plotting variables
Temp = zeros(len,N);                  %temp of the electrolyzer, [C]
U = zeros(len,1);                       %voltage/cell in each of the electrolyzer, [V]
I = zeros(len,1);                       %current in each electrolyzer, [A]
P = zeros(len,1);
I_den = zeros(len,1);                   %current density in the electrolyzer, [A/m^2]
nH2elout = zeros(len,1);                %hydrogen flowrate from each of the individual electrolyzer, [mol/s]
Telout = zeros(len,1);
Telin = zeros(len,1);
mBufferT = zeros(len,1);
Tw_out = zeros(len,1);
SpecEl = zeros(len,1);
PcompH2 = zeros(len,1);                 %compressor power for hydrogen, [watts]
PcompO2 = zeros(len,1);                 %compressor power for oxygen, [watts]
Qloss = zeros(len,1);                   %heat loss to surrounding in the electrolyzer, [watts]
Qgen = zeros(len,1);                    %heat generated in the electrolyzer, [watts]
Qlosslye = zeros(len,1);                %heat taken out by the lye from the electrolyzer, [watts]
P_net = zeros(len,1);                   %net power to the electrolyzer assembly, [watts]

%% Integrate plant over the time horizon

for i=1:len
    %i = timestamp
    %j = electrolyzer sequence
    
    r =  F('x0',x0,'z0',z0,'p',[V_El(i,:), qlye(i,:), q_cw(i,:), qH2O(i,:)]);
    x0 = full(r.xf);            %updating solution as new initial conditions
    z0 = full(r.zf);
    
    
    %% Storing values in plotting variables
    
    for j=1:N
        U(i,j) = full(r.zf(j));             %voltage/cell, [V]
        I(i,j) = full(r.zf(N+j));             %current, [A]
        P(i,j) = full(r.zf(2*N+j));             %power, [Watts]
        nH2elout(i,j) = full(r.zf(4*N+j));      %hydrogen production rate from individual electrolyzer, [mol/s]
        
        Temp(i,j)= full(r.xf(j));             %temperature of electrolyzers at all timestamps, [celsius]
        mBufferT(i,j) = full(r.xf(N+j));    %
        Tbtout(i,j) = full(r.xf(2*N+j));      %temperature of lye mixture at the exit of the buffer tank, [degC]
        Telin(i,j) = full(r.xf(3*N+j));       %temperature of lye going into the electrolyzer, [celsius]
        Tw_out(i,j) = full(r.xf(4*N+j));      %exit temperature of the cooling water, [celsius]
        
        
        I_den(i,j) = 0.1*I(i,j)/par.EL(j).A;    %current density, [mA/cm^2]
        
        
        V_H2(i,j) = nH2elout(i,j)*0.0224136*3600;%hydrogen production rate from individual electrolyzer, [Nm3/h]
        Ps(i,j) = P(i,j)/(1000*V_H2(i,j));%Specific electricity consumption, [kWh/Nm3]
    end
    
    P_net(i)=sum(P(i,:));
    
    if rem(i,100)==0
        disp(i)
    end
    
end

%% Plotting the results

figure()
subplot(2,1,1)
plot(V_H2)
xlabel('Time, s')
ylabel('H_2 production rate, [Nm^3/hr]')
legend('El 1','El 2', 'El 3')
grid on

subplot(2,1,2)
plot(Ps)
xlabel('Time, s')
ylabel('Specific electricity consumption, [kWh/Nm^3]')
ylim([4.3, 4.7])
legend('El 1','El 2', 'El 3')
grid on


figure()
subplot(2,1,1)
plot(Temp(:,1))
hold on
plot(Temp(:,2))
hold on
plot(Temp(:,3))
xlabel('Time, s')
ylabel('T_k, [ ^0C]')
ylim([69, 82])
grid on
subplot(2,1,2)
plot(Telout)
xlabel('Time, s')
ylabel('T_o_u_t,[ ^0 C]')
grid on

figure()
subplot(3,1,1)
plot(qlye(:,1))
xlabel('Time, s')
ylabel('Lye flowrate to El 1, g/s')
grid on
subplot(3,1,2)
plot(q_cw)
xlabel('Time, s')
ylabel('Cooling water flowrate, g/s')
grid on
subplot(3,1,3)
plot(Telin,'k')
hold on
plot(T_El_in_set,'r--')
xlabel('Time, s')
ylabel('T_E_l_ _i_n, [ ^0C]')
%ylim([63.5 65.5])
grid on

figure()
plot(mBufferT./1000)
ylabel('Mass of liquid in the buffer tank, [kg]')
xlabel('Time, s')
grid on
