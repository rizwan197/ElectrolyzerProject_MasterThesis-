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
num_hr = 1;                             %no. of hours
t0 = 1;                                 %start, [s)]
ts = 1;                                 %time step, [s]
tf = num_hr*60*60;                      %final, [s]
tsamp = t0:ts:tf;
len = length(tsamp);                    %number of simulation time steps
tstep = 200;
%% Initial guess for steady state optimization using IPOPT
%algebriac state variables('z')
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

%differential state variables('x')
T_k0 = 75*ones(1,par.N);
Psto_H20 = 25;      %initial H2 storage pressure (calculated from steady state solution) [bar]
Psto_O20 = 25;      %initial O2 storage pressure (calculated from steady state solution) [bar]
Mass_Bt0 = 1000000; %mass of the water in the buffer tank,[g]
T_El_in0 = 65;      %initial guess for the temperature of inlet lye into the electrolyzer, [deg C]
T_cw_out0 = 20;     %initial guess for the exit temperature of the cooling water leaving heat exchanger, [deg C]

z_guess = [u_k0 i_k0 P_k0 Feff_k0 nH2_k0 qH2Oloss_k0 nH2El_net0 nH2out_net0 nO2El_net0 nO2out_net0 T_El_out0];
x_guess = [T_k0 Psto_H20 Psto_O20 Mass_Bt0 T_El_in0 T_cw_out0];

%initial guess for input variables('u')
U_El_k_0 = 414.0301*ones(1,par.N);    %voltage across electrolyzers, [Volts]
q_lye_k_0 = 6648*ones(1,par.N);
q_cw_0 = 2.0698e4;                    %cooling water flow rate, [g/s]
zH2_0 = 0.4;
zO2_0 = 0.4;
q_H2O_0 = 324.2657;                    %total water lost during electrolysis, [grams/sec]

u_guess = [U_El_k_0 q_lye_k_0 q_cw_0 zH2_0 zO2_0 q_H2O_0]; 

%initial guess vector for the IPOPT
X_guess = [z_guess x_guess u_guess];


%% Solve the steady state optimization problem
[z0, x0, u0] = El_SteadyStateOptimization(N,X_guess);

T_El_in_set = x0(par.N+4);%setpoint for the temperature of lye entering the electrolyzer 
%Initial value of the MVs 
Vss = u0(1:par.N);
q_lyek = u0(par.N+1:2*par.N);
qf_cw = u0(2*par.N+1);
zH2 = u0(2*par.N+2);
zO2 = u0(2*par.N+3);
Qwater = u0(2*par.N+4);

%% Build the plant model
[xDiff, xAlg, input, eqnAlg, eqnDiff] = model(par.N);
x = [xAlg;xDiff];

%% Manipulated variables
%these are the degree of freedoms that we will utilise to control the system

V_El = zeros(len,N);                       %voltage across the electrolyzer, [Watt], len is the length of time vector
for j = 1:N
    V_El(1:end,j) = Vss(j)*1;               %incremental step change in common voltage across all electrolysers
    V_El(tstep:end,j)=Vss(j)*1;
end

qlye = zeros(len,N);                        %lye flowrate, [g/s]
for j = 1:N
    qlye(1:end,j) = q_lyek(j)*1;       %assumed same lye flowarate to all the electrolyzers
end
qlye(tstep:end,1) = q_lyek(2)*1.2;

q_cw = qf_cw*ones(len,1);                     %cooling water flow rate as a manipulated variable, [g/s]
q_cw(tstep:end) = qf_cw*1;                  %incremental step change in cooling water flowrate

ZH2 = zH2*ones(len,1);      %H2 valve displacement as a manipulated variable
%ZH2(tstep:end) = .2;                         %change in H2 valve displacement

ZO2 = zO2*ones(len,1);      %O2 valve displacement as a manipulated variable
%ZO2(tstep:end) = .7;                         %change in O2 valve displacement

qH2O = Qwater*ones(len,1);                  %flow rate of water added to buffer tank as a manipulated variable, [g/s]
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
Telin = zeros(len,1);
level = zeros(len,1);
SpecEl = zeros(len,1);
PcompH2 = zeros(len,1);                 %compressor power for hydrogen, [watts]
PcompO2 = zeros(len,1);                 %compressor power for oxygen, [watts]
Qloss = zeros(len,1);                   %heat loss to surrounding in the electrolyzer, [watts]
Qgen = zeros(len,1);                    %heat generated in the electrolyzer, [watts]
Qlosslye = zeros(len,1);                %heat taken out by the lye from the electrolyzer, [watts]
P_net = zeros(len,1);                   %net power to the electrolyzer assembly, [watts]

%% Build integrator to integrate plant over the time horizon

dae = struct('x',xDiff,'z',xAlg,'p',input,'ode',eqnDiff,'alg',eqnAlg);
F = integrator('F', 'idas', dae);

for i=1:len
    %i = timestamp
    %j = electrolyzer sequence
    
    r =  F('x0',x0,'z0',z0,'p',[V_El(i,:), qlye(i,:), q_cw(i), ZH2(i), ZO2(i), qH2O(i)]);
    x0 = full(r.xf);            %updating solution as new initial conditions
    z0 = full(r.zf);
    
    
    %% Storing values in plotting variables
    
    %calculation of compressor power
    PcompH2(i) = full(((r.zf(6*N+1)*par.Comp.k*par.Const.R*par.Comp.Tel)/(par.Comp.alpha*(par.Comp.k-1)))*(((r.xf(N+1)/par.Comp.Pel)^((par.Comp.k-1)/par.Comp.k))-1));
    PcompO2(i) = full(((r.zf(6*N+3)*par.Comp.k*par.Const.R*par.Comp.Tel)/(par.Comp.alpha*(par.Comp.k-1)))*(((r.xf(N+2)/par.Comp.Pel)^((par.Comp.k-1)/par.Comp.k))-1));
    %assuming same k and Tel for O2
    
    PstoH2(i) = full(r.xf(N+1));          %hydrogen storage pressure at all timestamps, [bar]
    PstoO2(i) = full(r.xf(N+2));          %oxygen storage pressure at all timestamps, [bar]
    nH2in(i) = full(r.zf(6*N+1));           %net hydrogen flow rate in to the storage at all timestamps, [mol/s]
    nH2out(i) = full(r.zf(6*N+2));          %net hydrogen flowrate out from the storage at all timestamps, [mol/s]
    nO2in(i) = full(r.zf(6*N+3));           %net oxygen flow rate in to the storage at all timestamps, [mol/s]
    nO2out(i) = full(r.zf(6*N+4));          %net oxygen flowrate out from the storage at all timestamps, [mol/s]
    Telout(i) = full(r.zf(6*N+5));            %temperature of lye after mixing before going to the buffer tank, [celsius]
    Telin(i) = full(r.xf(N+4));             %temperature of lye going into the electrolyzer, [celsius]
    Tw_out(i) = full(r.xf(N+5));          %exit temperature of the cooling water, [celsius]
    level(i) = full(r.xf(N+3));
    
    for j=1:N
        U(i,j) = full(r.zf(j));             %voltage/cell, [V]
        I(i,j) = full(r.zf(N+j));             %current, [A]
        P(i,j) = full(r.zf(2*N+j));             %power, [Watts]
        Temp(i,j)= full(r.xf(j));             %temperature of electrolyzers at all timestamps, [celsius]
        I_den(i,j) = 0.1*I(i,j)/par.EL(j).A;    %current density, [mA/cm^2]
        nH2elout(i,j) = full(r.zf(4*N+j));      %hydrogen production rate from individual electrolyzer, [mol/s]
        SpecEl(i,j) = (P(i,j)*(10^-6))./(nH2elout(i,j)*par.Const.MwtH2*(10^-6)*3600);   %specific electricity, [MWh/tonne H2]
        Qloss(i,j) = par.TherMo(j).A_surf*I(i,j)* par.EL(j).nc/1000*(par.TherMo(j).hc*(Temp(i,j)-par.EL(j).Ta) + ...
            par.sigma*par.em*((Temp(i,j)+273.15)^4-(par.EL(j).Ta+273.15)^4));             %heat loss to surrounding in the electrolyzer, [watts]
        Qgen(i,j) = par.EL(j).nc*(U(i,j)-par.EL(j).Utn)*I(i,j);                         %heat generated in the electrolyzer, [watts]
        Qlosslye(i,j) = qlye(i,j)*par.Const.CpLye*(Telin(i)-Temp(i,j));                 %heat taken out by the lye from the electrolyzer, [watts]
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
ylabel('H_2 production rate, Nm^3/hr')
legend('El 1','El 2', 'El 3')
grid on

subplot(2,1,2)
plot(Ps)
xlabel('Time, s')
ylabel('Specific electricity consumption, kWh/Nm^3')
legend('El 1','El 2', 'El 3')
grid on

figure()
subplot(2,1,1)
plot(PstoH2)
xlabel('Time, s')
ylabel('H_2 Storage pressure, bar')
grid on

figure()
subplot(2,1,1)
plot(Temp(:,1))
hold on
plot(Temp(:,2))
hold on
plot(Temp(:,3))
xlabel('Time, s')
ylabel('T_k, C')
grid on
subplot(2,1,2)
plot(Telout)
xlabel('Time, s')
ylabel('T_o_u_t, C')
grid on

figure()
subplot(3,1,1)
plot(qlye(:,1))
xlabel('Time, s')
ylabel('Lye flowrate to El 1')
grid on
subplot(3,1,2)
plot(q_cw)
xlabel('Time, s')
ylabel('Cooling duty')
grid on
subplot(3,1,3)
plot(Telin,'k')
hold on
plot(T_El_in_set,'r--')
xlabel('Time, s')
ylabel('T_E_l_ _i_n, C')
%ylim([63.5 65.5])
grid on