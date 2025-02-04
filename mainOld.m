%this is the old main file which was used for getting the responses for a
%step change in all the DOF. Use this with model.m
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
num_hr = .5;                           %no. of hours
t0 = 1;                                 %start, [s)]
ts = 1;                                 %time step, [s]
tf = num_hr*60*60;                      %final, [s]
tsamp = t0:ts:tf;
len = length(tsamp);                    %number of simulation time steps
tstep = 200;

%% Initial guess for steady state solution using IPOPT

%disturbance is total power
Pnet = 5.84e6;%1900e3*par.N; %6.3MW total input power

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

X_guess = [z_guess x_guess u_guess];

%% Solve the steady state problem
[z0, x0, u0, EXIT] = El_SteadyStateOptimization(N,X_guess,Pnet);

z_guess = z0;
x_guess = x0;
u_guess = u0;

% z0 = [Uk Ik Pk Feffk nH2k qH2Olossk nH2El_tot nH2out_tot nO2El_tot nO2out_tot T_el_out];
% x0 = [Tk PstoH2 PstoO2 massBt T_bt_out T_el_in T_CW_out];
% u0 = [Vss q_lyek qf_cw zH2 zO2 Qwater];

T_bt_out = x0(par.N+4);
T_El_in_set = x0(par.N+5);%setpoint for the temperature of lye entering the electrolyzer
T_cw_out = x0(par.N+6);

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
    Qgenk(nEl) = par.EL(nEl).nc*(z0(nEl)-par.EL(nEl).Utn)*z0(par.N+nEl)/1000;
    Qlossk(nEl) = par.TherMo(nEl).A_El*(par.TherMo(nEl).hc*(Tk(nEl)-par.EL(nEl).Ta) + par.sigma*par.em*((Tk(nEl)+273.15)^4-(par.EL(nEl).Ta+273.15)^4))/1000;
    Qlyelossk(nEl) = (q_lyek(nEl)*par.Const.CpLye*(T_El_in_set-Tk(nEl)))/1000;
    Qnet(nEl) = Qgenk(nEl)+Qlyelossk(nEl)-Qlossk(nEl);
end

%% Manipulated variables/Parameters for the dynamic simulation
%these are the degree of freedoms that we will utilise to control the system

V_El = zeros(len,N);              %voltage across the electrolyzer, [Watt], len is the length of time vector
for j = 1:N
    V_El(1:end,j) = Vss(j)*1;     %incremental step change in common voltage across all electrolysers
    %     V_El(tstep:end,j)=Vss(j)*1.02;
end

qlye = zeros(len,N);                   %lye flowrate, [g/s]
for j = 1:N
    qlye(1:end,j) = q_lyek(j)*1;       %assumed same lye flowarate to all the electrolyzers
end
% qlye(tstep:end,3) = q_lyek(3)*1.2;

q_cw = qf_cw*ones(len,1);                     %cooling water flow rate as a manipulated variable, [g/s]
% q_cw(tstep:end) = qf_cw*1.2;                  %incremental step change in cooling water flowrate

ZH2 = zH2*ones(len,1);      %H2 valve displacement as a manipulated variable
% ZH2(tstep:end) = .2;                         %change in H2 valve displacement

ZO2 = zO2*ones(len,1);      %O2 valve displacement as a manipulated variable
% ZO2(tstep:end) = .2;                         %change in O2 valve displacement

qH2O = Qwater*ones(len,1);  %flow rate of water added to buffer tank as a manipulated variable, [g/s]
qH2O(tstep:end)=Qwater*1.2;                  %incremental step change in the water flow rate

%% Initialize plotting variables
Temp = zeros(len,N);                  %temp of the electrolyzer, [C]
PstoH2 = zeros(len,1);                %H2 storage pressure, [bar]
PstoO2 = zeros(len,1);                %O2 storage pressure, [bar]
U = zeros(len,1);                       %voltage/cell in each of the electrolyzer, [V]
I = zeros(len,1);                       %current in each electrolyzer, [A]
P = zeros(len,1);
I_den = zeros(len,1);                   %current density in the electrolyzer, [A/m^2]
nH2in = zeros(len,1);                   %net hydrogen flow rate in to the storage, [mol/s]
nH2out = zeros(len,1);                  %net hydrogen flowrate out from the storage, [mol/s]
nH2elout = zeros(len,1);                %hydrogen flowrate from each of the individual electrolyzer, [mol/s]
nO2in = zeros(len,1);                   %net oxygen flow rate in to the storage, [mol/s]
nO2out = zeros(len,1);                  %net oxygen flowrate out from the storage, [mol/s]
Telout = zeros(len,1);
Tbtout = zeros(len,1);                  %temperature of lye mixture at the exit of the buffer tank, [degC]
Telin = zeros(len,1);
mBufferT = zeros(len,1);
Tw_out = zeros(len,1);
SpecEl = zeros(len,1);
PcompH2 = zeros(len,1);                 %compressor power for hydrogen, [watts]
PcompO2 = zeros(len,1);                 %compressor power for oxygen, [watts]
Qloss = zeros(len,N);                   %heat loss to surrounding in the electrolyzer, [watts]
Qgen = zeros(len,N);                    %heat generated in the electrolyzer, [watts]
Qlosslye = zeros(len,N);                %heat taken out by the lye from the electrolyzer, [watts]
P_net = zeros(len,1);                   %net power to the electrolyzer assembly, [watts]

 %% Build the plant model- without PI controller
[xDiff, xAlg, input, eqnAlg, eqnDiff, F] = model(par.N);

%% Integrate plant over the time horizon

for i=1:len
    %i = timestamp
    %j = electrolyzer sequence
    
    r =  F('x0',x0,'z0',z0,'p',[V_El(i,:), qlye(i,:), q_cw(i), ZH2(i), ZO2(i), qH2O(i)]);
    x0 = full(r.xf);            %updating solution as new initial conditions
    z0 = full(r.zf);
    
    
    %% Storing values in plotting variables
    
    %calculation of compressor power
    PcompH2(i) = full(((r.zf(6*N+1)*par.Comp.k*par.Const.R*par.Comp.Tel)/...
        (par.Comp.alpha*(par.Comp.k-1)))*(((r.xf(N+1)/par.Comp.Pel)^((par.Comp.k-1)/par.Comp.k))-1));
    PcompO2(i) = full(((r.zf(6*N+3)*par.Comp.k*par.Const.R*par.Comp.Tel)/...
        (par.Comp.alpha*(par.Comp.k-1)))*(((r.xf(N+2)/par.Comp.Pel)^((par.Comp.k-1)/par.Comp.k))-1));
    %assuming same k and Tel for O2
    
    
    nH2in(i) = full(r.zf(6*N+1));     %net hydrogen flow rate in to the storage at all timestamps, [mol/s]
    nH2out(i) = full(r.zf(6*N+2));    %net hydrogen flowrate out from the storage at all timestamps, [mol/s]
    nO2in(i) = full(r.zf(6*N+3));     %net oxygen flow rate in to the storage at all timestamps, [mol/s]
    nO2out(i) = full(r.zf(6*N+4));    %net oxygen flowrate out from the storage at all timestamps, [mol/s]
    Telout(i) = full(r.zf(6*N+5));    %temperature of lye after mixing before going to the buffer tank, [celsius]
    
    PstoH2(i) = full(r.xf(N+1));        %hydrogen storage pressure at all timestamps, [bar]
    PstoO2(i) = full(r.xf(N+2));        %oxygen storage pressure at all timestamps, [bar]
    mBufferT(i) = full(r.xf(N+3));      %mass of liquid in the buffer tank. [grams]
    Tbtout(i) = full(r.xf(N+4));        %temperature of lye mixture at the exit of the buffer tank, [degC]
    Telin(i) = full(r.xf(N+5));         %temperature of lye going into the electrolyzer, [celsius]
    Tw_out(i) = full(r.xf(N+6));        %exit temperature of the cooling water, [celsius]
    
    % z0 = [Uk Ik Pk Feffk nH2k qH2Olossk nH2El_tot nH2out_tot nO2El_tot nO2out_tot T_el_out];
    % x0 = [Tk PstoH2 PstoO2 massBt T_bt_out T_el_in T_CW_out];
    % u0 = [Vss q_lyek qf_cw zH2 zO2 Qwater];
    
    for j=1:N
        U(i,j) = full(r.zf(j));             %voltage/cell, [V]
        I(i,j) = full(r.zf(N+j));             %current, [A]
        P(i,j) = full(r.zf(2*N+j));             %power, [Watts]
        nH2elout(i,j) = full(r.zf(4*N+j));      %hydrogen production rate from individual electrolyzer, [mol/s]
        
        Temp(i,j)= full(r.xf(j));             %temperature of electrolyzers at all timestamps, [celsius]
        
        I_den(i,j) = 0.1*I(i,j)/par.EL(j).A;    %current density, [mA/cm^2]
        
        %     Qgenk(nEl) = par.EL(nEl).nc*(z0(nEl)-par.EL(nEl).Utn)*z0(par.N+nEl)/1000;
        %     Qlossk(nEl) = par.TherMo(nEl).A_El*(par.TherMo(nEl).hc*(Tk(nEl)-par.EL(nEl).Ta) + par.sigma*par.em*((Tk(nEl)+273.15)^4-(par.EL(nEl).Ta+273.15)^4))/1000;
        %     Qlyelossk(nEl) = (q_lyek(nEl)*par.Const.CpLye*(T_El_in_set-Tk(nEl)))/1000;
        %     Qnet(nEl) = Qgenk(nEl)+Qlyelossk(nEl)-Qlossk(nEl)
        Qloss(i,j) = par.TherMo(j).A_El*(par.TherMo(j).hc*(Temp(i,j)-par.EL(j).Ta) + ...
            par.sigma*par.em*((Temp(i,j)+273.15)^4-(par.EL(j).Ta+273.15)^4));
        %heat loss to surrounding in the electrolyzer, [watts]
        
        Qgen(i,j) = par.EL(j).nc*(U(i,j)-par.EL(j).Utn)*I(i,j);
        %heat generated in the electrolyzer, [watts]
        
        Qlosslye(i,j) = qlye(i,j)*par.Const.CpLye*(Telin(i)-Temp(i,j));
        %heat taken out by the lye from the electrolyzer, [watts]
        Qnet(i,j) = Qgen(i,j)+Qlosslye(i,j)-Qloss(i,j);
        
        V_H2(i,j) = nH2elout(i,j)*0.0224136*3600;%hydrogen production rate from individual electrolyzer, [Nm3/h]
        Ps(i,j) = P(i,j)/(1000*V_H2(i,j));%Specific electricity consumption, [kWh/Nm3]
    end
    
    

    P_net(i)=sum(P(i,:));
    
    if rem(i,1000)==0
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
subplot(6,2,1)
plot(qlye)
xlabel('Time, s')
ylabel('Lye flowrate, g/s')
subplot(6,2,3)
plot(V_El)
xlabel('Time, s')
ylabel('Cell voltage')
subplot(6,2,5)
plot(q_cw./1000)
xlabel('Time, s')
ylabel('q_{cw}, kg/s')
subplot(6,2,7)
plot(ZH2)
xlabel('Time, s')
ylabel('z_{H_2}')
subplot(6,2,9)
plot(ZO2)
xlabel('Time, s')
ylabel('z_{O_2}')
subplot(6,2,11)
plot(qH2O./1000)
xlabel('Time, s')
ylabel('q_{H_2O}, kg/s')
subplot(6,2,2)
plot(Temp)
xlabel('Time, s')
ylabel('T_k, [ ^0C]')
subplot(6,2,4)
plot(I_den)
xlabel('Time, s')
ylabel('I_{den}, mA/cm^2')
subplot(6,2,6)
plot(Telin)
xlabel('Time, s')
ylabel('T_{El,in}, [ ^0C]')
subplot(6,2,8)
plot(PstoH2)
xlabel('Time, s')
ylabel('p_{sto} H_2, bar')
subplot(6,2,10)
plot(PstoO2)
xlabel('Time, s')
ylabel('p_{sto} O_2, bar')
subplot(6,2,12)
plot(mBufferT./1000)
ylabel('Mass_{bt}, [kg]')
xlabel('Time, s')


figure()
subplot(2,1,1)
plot(qH2O./1000)
xlabel('Time, s')
ylabel('q_{H_2O}, kg/s')
subplot(2,1,2)
plot(mBufferT./1000)
ylabel('Mass_{bt}, [kg]')
xlabel('Time, s')


%% Plot for checking Qnet in each electrolyzer
% figure()
% subplot(2,3,1)
% plot(Qnet(:,1))
% xlabel('Time, s')
% ylabel('Q_{net,1}')
% subplot(2,3,2)
% plot(Qnet(:,2))
% xlabel('Time, s')
% ylabel('Q_{net,2}')
% subplot(2,3,3)
% plot(Qnet(:,3))
% xlabel('Time, s')
% ylabel('Q_{net,3}')
% subplot(2,3,4)
% plot(Temp(:,1))
% xlabel('Time, s')
% ylabel('T_{el,1}')
% subplot(2,3,5)
% plot(Temp(:,2))
% xlabel('Time, s')
% ylabel('T_{el,2}')
% subplot(2,3,6)
% plot(Temp(:,3))
% xlabel('Time, s')
% ylabel('T_{el,3}')
%% Creating the data file
% save('data_Step_qlye2_40hr')
% load data_q1step3000_MVQcool_12hr
