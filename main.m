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

%Disturbance is total power
Psupp = [5.8e6];

%% Initialize plotting variables
Temp = zeros(len*length(Psupp),N);                  %temp of the electrolyzer, [C]
PstoH2 = zeros(len*length(Psupp),1);                %H2 storage pressure, [bar]
PstoO2 = zeros(len*length(Psupp),1);                %O2 storage pressure, [bar]
U = zeros(len*length(Psupp),1);                       %voltage/cell in each of the electrolyzer, [V]
I = zeros(len*length(Psupp),1);                       %current in each electrolyzer, [A]
P = zeros(len*length(Psupp),1);
I_den = zeros(len*length(Psupp),1);                   %current density in the electrolyzer, [A/m^2]
nH2in = zeros(len*length(Psupp),1);                   %net hydrogen flow rate in to the storage, [mol/s]
nH2out = zeros(len*length(Psupp),1);                  %net hydrogen flowrate out from the storage, [mol/s]
nH2elout = zeros(len*length(Psupp),1);                %hydrogen flowrate from each of the individual electrolyzer, [mol/s]
nO2in = zeros(len*length(Psupp),1);                   %net oxygen flow rate in to the storage, [mol/s]
nO2out = zeros(len*length(Psupp),1);                  %net oxygen flowrate out from the storage, [mol/s]
Telout = zeros(len*length(Psupp),1);
Tbtout = zeros(len*length(Psupp),1);                  %temperature of lye mixture at the exit of the buffer tank, [degC]
Telin = zeros(len*length(Psupp),1);
mBufferT = zeros(len*length(Psupp),1);
Tw_out = zeros(len*length(Psupp),1);
SpecEl = zeros(len*length(Psupp),1);
PcompH2 = zeros(len*length(Psupp),1);                 %compressor power for hydrogen, [watts]
PcompO2 = zeros(len*length(Psupp),1);                 %compressor power for oxygen, [watts]
Qloss = zeros(len*length(Psupp),1);                   %heat loss to surrounding in the electrolyzer, [watts]
Qgen = zeros(len*length(Psupp),1);                    %heat generated in the electrolyzer, [watts]
Qlosslye = zeros(len*length(Psupp),1);                %heat taken out by the lye from the electrolyzer, [watts]
P_net = zeros(len*length(Psupp),1);                   %net power to the electrolyzer assembly, [watts]

%MVs for the control of the process
V_El = zeros(len*length(Psupp),1);                  %voltage across the electrolyzer, [Watt], len is the length of time vector
qlye = zeros(len*length(Psupp),N);                  %lye flowrate, [g/s]
q_cw = zeros(len*length(Psupp),1);                  %cooling water flowrate
%%
for distInt = 1:length(Psupp)
Pnet = Psupp(distInt);%total input power

%% Initial guess for steady state solution using IPOPT

[z_guess,x_guess,u_guess] = init0(par.N,Pnet);


%% Solve the steady state problem
X_guess = [z_guess x_guess u_guess];

z0 = z_guess;
x0 = x_guess;
u0 = u_guess;

%Initial value of the MVs
Vss = u0(1);%electrolyzer voltage
q_lyek = u0(2:par.N+1);
qlye_kgs = q_lyek/1000;
qf_cw = u0(par.N+2);
qcw_kgs = qf_cw/1000;
zH2 = u0(par.N+3);
zO2 = u0(par.N+4);
Qwater = u0(par.N+5);

nH2 = sum(z0(4*par.N+1:5*par.N));
Tk = x0(1:par.N);
T_El_in = x0(par.N+5);      %temperature of lye entering the electrolyzer
T_cw_out = x0(par.N+6);
T_bt_out = x0(par.N+4);

Pcons = sum(z0(2*par.N+1:3*par.N));

Iden = 0.1*z0(par.N+1:2*par.N)./par.EL(1).A;
V_H2_ini = z0(4*par.N+1:5*par.N)*0.0224136*3600;
% for nEl = 1:par.N
%     Qgenk(nEl) = par.EL(nEl).nc*(z0(nEl)-par.EL(nEl).Utn)*z0(par.N+nEl)/1000;
%     Qlossk(nEl) = par.TherMo(nEl).A_El*(par.TherMo(nEl).hc*(Tk(nEl)-par.EL(nEl).Ta) + par.sigma*par.em*((Tk(nEl)+273.15)^4-(par.EL(nEl).Ta+273.15)^4))/1000;
% end


%% Manipulated variables/Parameters for the dynamic simulation
%these are the degree of freedoms that we will utilise to control the system


for i=1+len*(distInt-1):len*distInt

for j = 1:N
    qlye(i,j) = q_lyek(j)*1;       %assumed same lye flowrate to all the electrolyzers
end
% qlye(tstep+500:end,1) = q_lyek(1)*1.2;

q_cw(i) = qf_cw;                     %cooling water flow rate as a manipulated variable, [g/s]
% q_cw(tstep+500:end) = qf_cw*1.2;                  %incremental step change in cooling water flowrate
end


%% Build the plant model- with PI controller

%PI controller for states with zero SS effect
MV_0 = [zH2 zO2 Qwater];
tauC = 50;
K = [-1e-3, -1e-3, 0.975];
Kc = 1./(K*tauC);
tauI = 1*4*tauC*ones(1,3);

%proportional controller for total power
PC.u0 = Vss;
PC.Kc = 8.5/622000;

[xDiff, xAlg, input, eqnAlg, eqnDiff, F] = modelPI(par.N,MV_0,Kc,tauI,PC);
x0 = [x0,0,0,0];%input vector for the differential states including the eint terms

%define the setpoint trajectory
pstoH2set = x0(par.N+1)*ones(len*length(Psupp),1);
% pstoH2set = 24.77*ones(len*length(Psupp),1);
% pstoH2set(tstep:end) = 24.75;
pstoO2set = x0(par.N+2)*ones(len*length(Psupp),1);
% pstoO2set = 24.77*ones(len*length(Psupp),1);
% pstoO2set(tstep:end) = 24.75;
Mass_Btset = x0(par.N+3)*ones(len*length(Psupp),1); %setpoint for the mass in the buffer tank
% Mass_Btset = 3000*ones(len*length(Psupp),1);
% Mass_Btset(tstep:end) = 3500000;

Pnet_set = Pnet;%setpoint for the power

%% Simulate the plant

for i=1+len*(distInt-1):len*distInt
    %i = timestamp
    %j = electrolyzer sequence
    
    r =  F('x0',x0,'z0',z0,'p',[Pnet_set, qlye(i,:), q_cw(i), pstoH2set(i), pstoO2set(i), Mass_Btset(i)]);
    x0 = full(r.xf);            %updating solution as new initial conditions
    z0 = full(r.zf);
    
    
    %% Storing values in plotting variables
    
    nH2in(i) = full(r.zf(6*N+1));     %net hydrogen flow rate in to the storage at all timestamps, [mol/s]
    nH2out(i) = full(r.zf(6*N+2));    %net hydrogen flowrate out from the storage at all timestamps, [mol/s]
    nO2in(i) = full(r.zf(6*N+3));     %net oxygen flow rate in to the storage at all timestamps, [mol/s]
    nO2out(i) = full(r.zf(6*N+4));    %net oxygen flowrate out from the storage at all timestamps, [mol/s]
    Telout(i) = full(r.zf(6*N+5));    %temperature of lye after mixing before going to the buffer tank, [celsius]
    
    PstoH2(i) = full(r.xf(N+1));        %hydrogen storage pressure at all timestamps, [bar]
    PstoO2(i) = full(r.xf(N+2));        %oxygen storage  pressure at all timestamps, [bar]
    mBufferT(i) = full(r.xf(N+3));      %mass of liquid in the buffer tank. [grams]
    Tbtout(i) = full(r.xf(N+4));        %temperature of lye mixture at the exit of the buffer tank, [degC]
    Telin(i) = full(r.xf(N+5));         %temperature of lye going into the electrolyzer, [celsius]
    Tw_out(i) = full(r.xf(N+6));        %exit temperature of the cooling water, [celsius]
    eint_pstoH2(i) = full(r.xf(N+7));   %integrated error term for pstoH2
    eint_pstoO2(i) = full(r.xf(N+8));   %integrated error term for pstoO2
    eint_Mbt(i) = full(r.xf(N+9));      %integrated error term for the mass of the buffer tank
    
    for j=1:N
        U(i,j) = full(r.zf(j));             %voltage/cell, [V]
        I(i,j) = full(r.zf(N+j));             %current, [A]
        P(i,j) = full(r.zf(2*N+j));             %power, [Watts]
        nH2elout(i,j) = full(r.zf(4*N+j));      %hydrogen production rate from individual electrolyzer, [mol/s]
        
        Temp(i,j)= full(r.xf(j));             %temperature of electrolyzers at all timestamps, [celsius]
        
        I_den(i,j) = 0.1*I(i,j)/par.EL(j).A;    %current density, [mA/cm^2]
        
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
    
    %% Calculate the input trajectory with PI controller
    % for storage pressure of the hydrogen tank
    Kc_pstoH2PI = Kc(1);                    %controller gain
    taui_pstoH2PI = tauI(1);                %integral time constant
    e_pstoH2 = (PstoH2(i) - pstoH2set(i));  %error term
    ZH2(i) = PIcontroller(zH2,Kc_pstoH2PI,taui_pstoH2PI,e_pstoH2,eint_pstoH2(i));
    
    
    %for storage pressure of the oxygen tank
    Kc_pstoO2PI = Kc(2);                    %controller gain
    taui_pstoO2PI = tauI(2);                %integral time constant
    e_pstoO2 = (PstoO2(i) - pstoO2set(i));
    ZO2(i) = PIcontroller(zO2,Kc_pstoO2PI,taui_pstoO2PI,e_pstoO2,eint_pstoO2(i));

    
    %for mass in the buffer tank
    Kc_MassBtPI = Kc(3);
    taui_MassBtPI = tauI(3);
    e_Mbt = mBufferT(i) - Mass_Btset(i);
    q_H2O(i) = PIcontroller(Qwater,Kc_MassBtPI,taui_MassBtPI,e_Mbt,eint_Mbt(i));
    %%
    P_net(i)=sum(P(i,:)); 
    
    %controller for MV-CV pairing Vel-Pnet
    V_El(i) = PC.u0 - PC.Kc*(P_net(i)-Pnet_set);
    
    if rem(i,1000)==0
        disp(i)
    end
    
end

end

%% Plotting the results

% figure()
% subplot(2,1,1)
% plot(V_H2)
% xlabel('Time, s')
% ylabel('H_2 production rate, [Nm^3/hr]')
% legend('El 1','El 2', 'El 3')
% grid on
% subplot(2,1,2)
% plot(Ps)
% xlabel('Time, s')
% ylabel('Specific electricity consumption, [kWh/Nm^3]')
% ylim([4.3, 4.7])
% legend('El 1','El 2', 'El 3')
% grid on

figure()
subplot(6,2,1)
plot(qlye)
xlabel('Time, s')
ylabel('Lye flowrate, g/s')
subplot(6,2,3)
plot(U)
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
plot(q_H2O./1000)
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

figure()
subplot(2,1,1)
plot(V_El)
xlabel('Time, s')
ylabel('Cell voltage')
subplot(2,1,2)
plot(I_den)
xlabel('Time, s')
ylabel('I_{den}, mA/cm^2')

figure()
plot(P_net./1e6)
xlabel('Time, s')
ylabel('P_{net},MW')

%% Creating the data file
% save('data_closeloop_Velstep_40hr')
% load data_q1step3000_MVQcool_12hr
