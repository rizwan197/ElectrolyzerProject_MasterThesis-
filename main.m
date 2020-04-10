clc
clear
close all

%Solving the dynamic optimization problem for Flowsheet 1 State 2
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
Tp = 900; %prediction horizon for MPC
Tc = 900; %control horizon for MPC
tstep = 200;

%% Initial guess for steady state RTO

%disturbance is total power
Pnet = 5.6e6;%1900e3*par.N; %6.3MW total input power

X_guess = initializeRTO(N,Pnet);

%% Solve the steady state RTO, economic objective function
[z0, x0, u0] = El_RTO(N,X_guess,Pnet);
% T_El_in_set = x0(par.N+5);%setpoint for the temperature of lye entering the electrolyzer 

%Initial value of the MVs 
Vss = u0(1:par.N);
q_lyek = u0(par.N+1:2*par.N);
qlye_kgs = q_lyek/1000; 
qf_cw = u0(2*par.N+1);
qcw_kgs = qf_cw/1000;
zH2 = u0(2*par.N+2);
zO2 = u0(2*par.N+3);
Qwater = u0(2*par.N+4);

Pcons = sum(z0(2*par.N+1:3*par.N));
Iden = 0.1*z0(par.N+1:2*par.N)./par.EL(1).A;
Tk = x0(1:par.N);
V_H2_ini = z0(4*par.N+1:5*par.N)*0.0224136*3600;
%Check for the RTO solution
% row = [Pcons/1e6,qlye_kgs,qcw_kgs,Iden,Tk,T_El_in_set,V_H2_ini];

%% Build the plant model
[xDiff, xAlg, input, eqnAlg, eqnDiff, F] = model(par.N);

%Add OCP for the MPC here 



%% Manipulated variables
%these are the degree of freedoms that we will utilise to control the system

V_El = zeros(len,N);              %voltage across the electrolyzer, [Watt], len is the length of time vector
for j = 1:N
    V_El(1:end,j) = Vss(j)*1;     %incremental step change in common voltage across all electrolysers
%     V_El(tstep:end,j)=Vss(j)*1;
end

qlye = zeros(len,N);                   %lye flowrate, [g/s]
for j = 1:N
    qlye(1:end,j) = q_lyek(j)*1;       %assumed same lye flowarate to all the electrolyzers
end
% qlye(tstep:end,2) = q_lyek(2)*1.2;

q_cw = qf_cw*ones(len,1);                     %cooling water flow rate as a manipulated variable, [g/s]
% q_cw(tstep:end) = qf_cw*1.1;                  %incremental step change in cooling water flowrate

ZH2 = zH2*ones(len,1);      %H2 valve displacement as a manipulated variable
%ZH2(tstep:end) = .2;                         %change in H2 valve displacement

ZO2 = zO2*ones(len,1);      %O2 valve displacement as a manipulated variable
%ZO2(tstep:end) = .7;                         %change in O2 valve displacement

qH2O = Qwater*ones(len,1);  %flow rate of water added to buffer tank as a manipulated variable, [g/s]
%qH2O(tstep:end)=Qwater*1.2;                  %incremental step change in the water flow rate

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
Qloss = zeros(len,1);                   %heat loss to surrounding in the electrolyzer, [watts]
Qgen = zeros(len,1);                    %heat generated in the electrolyzer, [watts]
Qlosslye = zeros(len,1);                %heat taken out by the lye from the electrolyzer, [watts]
P_net = zeros(len,1);                   %net power to the electrolyzer assembly, [watts]

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
    
    PstoH2(i) = full(r.xf(N+1));      %hydrogen storage pressure at all timestamps, [bar]
    PstoO2(i) = full(r.xf(N+2));      %oxygen storage pressure at all timestamps, [bar]
    mBufferT(i) = full(r.xf(N+3));    %
    Tbtout(i) = full(r.xf(N+4));      %temperature of lye mixture at the exit of the buffer tank, [degC]
    Telin(i) = full(r.xf(N+5));       %temperature of lye going into the electrolyzer, [celsius]
    Tw_out(i) = full(r.xf(N+6));      %exit temperature of the cooling water, [celsius]
    
    
    for j=1:N
        U(i,j) = full(r.zf(j));             %voltage/cell, [V]
        I(i,j) = full(r.zf(N+j));             %current, [A]
        P(i,j) = full(r.zf(2*N+j));             %power, [Watts]
        nH2elout(i,j) = full(r.zf(4*N+j));      %hydrogen production rate from individual electrolyzer, [mol/s]
        
        Temp(i,j)= full(r.xf(j));             %temperature of electrolyzers at all timestamps, [celsius]
        
        I_den(i,j) = 0.1*I(i,j)/par.EL(j).A;    %current density, [mA/cm^2]
        
        Qloss(i,j) = par.TherMo(j).A_El*(par.TherMo(j).hc*(Temp(j)-par.EL(j).Ta) + ...
            par.sigma*par.em*((Temp(j)+273.15)^4-(par.EL(j).Ta+273.15)^4)); 
        %heat loss to surrounding in the electrolyzer, [watts]
        
        Qgen(i,j) = par.EL(j).nc*(U(i,j)-par.EL(j).Utn)*I(i,j);             
        %heat generated in the electrolyzer, [watts]
        
        Qlosslye(i,j) = qlye(i,j)*par.Const.CpLye*(Telin(i)-Temp(i,j));                 
        %heat taken out by the lye from the electrolyzer, [watts]
        
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
plot(PstoH2)
xlabel('Time, s')
ylabel('H_2 Storage pressure, [bar]')
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
% hold on
% plot(T_El_in_set,'r--')
xlabel('Time, s')
ylabel('T_E_l_ _i_n, [ ^0C]')
%ylim([63.5 65.5])
grid on

figure()
plot(mBufferT./1000)
ylabel('Mass of liquid in the buffer tank, [kg]')
xlabel('Time, s')
grid on

%% Creating the data file
% save('data_Step_qlye2_40hr')
% load data_q1step3000_MVQcool_12hr
