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
num_hr = 3;                           %no. of hours
t0 = 1;                                 %start, [s)]
ts = 1;                                 %time step, [s]
tf = num_hr*60*60;                      %final, [s]
tsamp = t0:ts:tf;
len = length(tsamp);                    %number of simulation time steps
tstep = 200;

%% Initial guess for steady state solution using IPOPT

%disturbance is total power
P_inp = 5.5e6; % total input power

[z_guess,x_guess,u_guess] = init0(N,P_inp);
X_guess = [z_guess x_guess u_guess];

%% Solve the steady state problem
z0 = z_guess;
x0 = x_guess;
u0 = u_guess;

T_El_in = x0(par.N+5);%setpoint for the temperature of lye entering the electrolyzer
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

Iden = 0.1*z0(par.N+1:2*par.N)./par.EL(1).A;
V_H2_ini = z0(4*par.N+1:5*par.N)*0.0224136*3600;

%SOC for the qcw for 3.4<=P_inp<5.1MW c=Hymeas, here H = [-0.9726 0.2326],
%measurements = [Telin_set T_cwout_set]; these need to be updated online from RTO above
SOC.H = [-0.9726 0.2326];
SOC.ymeas_set = [67.57 72.94]';%at P_inp = 5MW, requires online updation

SOC.c_set = SOC.H*SOC.ymeas_set;
SOC.ymeas0 = [T_El_in T_cw_out]';

%initial value of the CVs
SOC.c(1) = SOC.H*SOC.ymeas0;
Tk(1,:) = x0(1:par.N);
Iden1(1) = Iden(1);
Pcons(1) = sum(z0(2*par.N+1:3*par.N));
T_ElinC(1) = T_El_in;
% T_ElinC(1) = max(Tk(1,:))-T_El_in;

% for nEl = 1:par.N
%     Qgenk(nEl) = par.EL(nEl).nc*(z0(nEl)-par.EL(nEl).Utn)*z0(par.N+nEl)/1000;
%     Qlossk(nEl) = par.TherMo(nEl).A_El*(par.TherMo(nEl).hc*(Tk(nEl)-par.EL(nEl).Ta) + par.sigma*par.em*((Tk(nEl)+273.15)^4-(par.EL(nEl).Ta+273.15)^4))/1000;
% end

%% Manipulated variables/Parameters for the dynamic simulation
%these are the degree of freedoms that we will utilise to control the system

% V_El = zeros(len,N);              %voltage across the electrolyzer, [Watt], len is the length of time vector
% for j = 1:N
%     V_El(1:end,j) = Vss(j)*1;     %incremental step change in common voltage across all electrolysers
% %         V_El(tstep+500:end,j)=Vss(j)*1.02;
% end
%
% qlye = zeros(len,N);                   %lye flowrate, [g/s]
% for j = 1:N
%     qlye(1:end,1) = 10e3;%q_lyek(1)*1;       %assumed same lye flowarate to all the electrolyzers
% end
% qlye(tstep+500:end,1) = q_lyek(1)*1.2;

% q_cw = qf_cw*ones(len,1);                     %cooling water flow rate as a manipulated variable, [g/s]
% q_cw(tstep+500:end) = qf_cw*1.2;                  %incremental step change in cooling water flowrate

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
P_net = 4.5e6*ones(len,1);              %net power consumed by the electrolyzer assembly, [watts]

%% Build the plant model with PI controller (Sigurd's book method)
MV_0 = [zH2 zO2 Qwater];
Kc = 1.*[-20 -20 0.021];
tauI = 1.*[200 200 200];
[xDiff, xAlg, input, eqnAlg, eqnDiff, F] = modelPI(par.N,MV_0,Kc,tauI);
x0 = [x0,0,0,0];%input vector for the differential states including the eint terms
%define the setpoint trajectory
pstoH2set = x0(par.N+1)*ones(len,1);
% pstoH2set(tstep:end) = 24.75;
pstoO2set = x0(par.N+2)*ones(len,1);
% pstoO2set(tstep:end) = 24.75;
Mass_Btset = x0(par.N+3)*ones(len,1); %setpoint for the mass in the buffer tank
% Mass_Btset(tstep:end) = 3500000;
%% PI controller for rest of the states
%% Region 1; P_inp>=5.84; no unconstrained DOF

% %fix qlye1 = 10kg/s
% qlye = zeros(len,N);
% qlye(1:end,1) = 10e3;%q_lyek(1)*1;
% 
% % %for pairing qcw-Iden1
% % Iden1C.u0 = qf_cw;
% % Iden1C.tauC = 100;
% % Iden1C.k = -1.465e-4;
% % Iden1C.tau1 = 15300;
% % Iden1C.Kc = (1/Iden1C.k)*(Iden1C.tau1/Iden1C.tauC);
% % Iden1C.tauI = min(Iden1C.tau1,4*Iden1C.tauC);
% % Iden1C.Ki = Iden1C.Kc/Iden1C.tauI;
% % Iden1C.set = 198.5;
% % Iden1C.err0 = 0;
% 
% %for pairing Vel-Iden1
% Iden1C.u0 = max(Vss);
% Iden1C.Kc = 8.5/16.5;
% Iden1C.set = 198.5;
% 
% %for pairing qlye2-T2
% T2C.u0 = q_lyek(2);
% T2C.tauC = 100;
% T2C.k = -9.712e-3;
% T2C.tau1 = 7800;
% T2C.Kc = (1/T2C.k)*(T2C.tau1/T2C.tauC);
% T2C.tauI = min(T2C.tau1,4*T2C.tauC);
% T2C.Ki = T2C.Kc/T2C.tauI;
% T2C.set = 80;
% T2C.err0 = 0;
% 
% %for pairing qlye3-T3
% T3C.u0 = q_lyek(3);
% T3C.tauC = 100;
% T3C.k = -0.01166;
% T3C.tau1 = 8900;
% T3C.Kc = (1/T3C.k)*(T3C.tau1/T3C.tauC);
% T3C.tauI = min(T3C.tau1,4*T3C.tauC);
% T3C.Ki = T3C.Kc/T3C.tauI;
% T3C.set = 80;
% T3C.err0 = 0;
% 
% % %pairing for Vel-(Telin)
% % TElinC.u0 = max(Vss);
% % TElinC.tauC = 100;
% % TElinC.kDash = 9.05e-5;
% % TElinC.Kc = (1/TElinC.kDash)*(1/TElinC.tauC);
% % TElinC.tauI = 4*TElinC.tauC;
% % TElinC.Ki = TElinC.Kc/TElinC.tauI;
% % TElinC.set = 50;
% % TElinC.err0 = 0;
% 
% %for pairing qcw-Telin
% TElinC.u0 = qf_cw;
% TElinC.Kc = -63000e-4;
% TElinC.set = 50;

%% Region 2; 5.1MW<=P_inp<5.84

%for pairing qlye1-Iden1
Iden1C.u0 = q_lyek(1);
Iden1C.tauC = 100;
Iden1C.k = -1.1e-3;
Iden1C.tau1 = 1600;
Iden1C.Kc = (1/Iden1C.k)*(Iden1C.tau1/Iden1C.tauC);
Iden1C.tauI = min(Iden1C.tau1,4*Iden1C.tauC);
Iden1C.Ki = Iden1C.Kc/Iden1C.tauI;
Iden1C.set = 198.5;
Iden1C.err0 = 0;

%for pairing qlye2-T2
T2C.u0 = q_lyek(2);
T2C.tauC = 100;
T2C.k = -9.712e-3;
T2C.tau1 = 7800;
T2C.Kc = (1/T2C.k)*(T2C.tau1/T2C.tauC);
T2C.tauI = min(T2C.tau1,4*T2C.tauC);
T2C.Ki = T2C.Kc/T2C.tauI;
T2C.set = 80;
T2C.err0 = 0;

%for pairing qlye3-T3
T3C.u0 = q_lyek(3);
T3C.tauC = 100;
T3C.k = -0.01166;
T3C.tau1 = 8900;
T3C.Kc = (1/T3C.k)*(T3C.tau1/T3C.tauC);
T3C.tauI = min(T3C.tau1,4*T3C.tauC);
T3C.Ki = T3C.Kc/T3C.tauI;
T3C.set = 80;
T3C.err0 = 0;

%for pairing Vel-Pnet
PC.u0 = max(Vss);
PC.Kc = 8.5/622000;
PC.set = P_inp;

%for pairing qcw-TElin; P controller
TElinC.u0 = qf_cw;
TElinC.Kc = -63000e-4;
TElinC.set = 56.69;%from RTO layer above

% %for pairing unconstrained qcw-SOC
% SOC.u0 = qf_cw;
% SOC.tauC = 100;
% SOC.k = 0.034;
% SOC.tau1 = 500;%2750;
% SOC.theta = 50;
% SOC.Kc = (1/SOC.k)*(SOC.tau1/(SOC.tauC+SOC.theta));
% SOC.tauI = min(SOC.tau1,4*(SOC.tauC+SOC.theta));
% SOC.Ki = SOC.Kc/SOC.tauI;
% SOC.set = SOC.c_set;
% SOC.err0 = 0;

%% Region 3; 3.4<=P_inp<5.1MW
% %for pairing qlye1-T1
% T1C.u0 = q_lyek(1);
% T1C.tauC = 100;
% T1C.k = -1.358e-3;
% T1C.tau1 = 2900;
% T1C.Kc = (1/T1C.k)*(T1C.tau1/T1C.tauC);
% T1C.tauI = min(T1C.tau1,4*T1C.tauC);
% T1C.Ki = T1C.Kc/T1C.tauI;
% T1C.set = 80;
% T1C.err0 = 0;
% 
% %for pairing qlye2-T2
% T2C.u0 = q_lyek(2);
% T2C.tauC = 100;
% T2C.k = -9.712e-3;
% T2C.tau1 = 7800;
% T2C.Kc = (1/T2C.k)*(T2C.tau1/T2C.tauC);
% T2C.tauI = min(T2C.tau1,4*T2C.tauC);
% T2C.Ki = T2C.Kc/T2C.tauI;
% T2C.set = 80;
% T2C.err0 = 0;
% 
% %for pairing qlye3-T3
% T3C.u0 = q_lyek(3);
% T3C.tauC = 100;
% T3C.k = -0.01166;
% T3C.tau1 = 8900;
% T3C.Kc = (1/T3C.k)*(T3C.tau1/T3C.tauC);
% T3C.tauI = min(T3C.tau1,4*T3C.tauC);
% T3C.Ki = T3C.Kc/T3C.tauI;
% T3C.set = 80;
% T3C.err0 = 0;
% 
% %for pairing Vel-Pnet
% PC.u0 = max(Vss);
% PC.Kc = 8.5/622000;
% PC.set = P_inp;
% 
% %for pairing qcw-Telin
% %PI controller
% % TElinC.tauC = 100;
% % TElinC.tau1 = 14300;
% % TElinC.Kc = (1/TelinC.k)*(TelinC.tau1/TelinC.tauC);
% % TElinC.tauI = min(TelinC.tau1,4*TelinC.tauC);
% % TElinC.Ki = TelinC.Kc/TelinC.tauI;
% % TElinC.set = 66.28;%from RTO layer above
% % TElinC.err0 = 0;
% %P controller
% TElinC.u0 = qf_cw;
% TElinC.Kc = -63000;
% TElinC.set = 66.28;%from RTO layer above

% %for pairing uncontrained DOF; qcw-SOC (T_elin,Tcw_out)
% SOC.u0 = qf_cw;
% SOC.tauC = 100;
% SOC.k = 0.034;
% SOC.tau1 = 500;%2750;
% SOC.theta = 50;
% SOC.Kc = (1/SOC.k)*(SOC.tau1/(SOC.tauC+SOC.theta));
% SOC.tauI = min(SOC.tau1,4*(SOC.tauC+SOC.theta));
% SOC.Ki = SOC.Kc/SOC.tauI;
% SOC.set = SOC.c_set;
% SOC.err0 = 0;

%% Integrate plant over the time horizon
for i=1:len
    
    %% Region 1; P_inp>5.84
    
%     Iden1C.err = Iden1C.set - Iden1(i);
%     T2C.err = T2C.set - Tk(i,2);
%     T3C.err = T3C.set - Tk(i,3);
%     TElinC.err = TElinC.set - T_ElinC(i);
%     
% %             %PI controller qcw-Iden1
% %                 q_cw(i) = min(80e3,max(1e3,Iden1C.u0 + Iden1C.Kc*(Iden1C.err-Iden1C.err0) + Iden1C.Ki*Iden1C.err*ts));%0.01<=qcw<=80
% %                 Iden1C.err0 = Iden1C.err;
% %                 Iden1C.u0 = q_cw(i);
%     
%     %Proportional controller Vel-Iden1
%     VEl(i) = Iden1C.u0 + Iden1C.Kc*(Iden1C.err);
%     Iden1C.u0 = VEl(i);
%     V_El(i,:) = Iden1C.u0*ones(1,par.N);
%     
%     %PI controller qlye2-T2
%     qlye(i,2) = min(10e3,max(0.5e3,T2C.u0 + T2C.Kc*(T2C.err-T2C.err0) + T2C.Ki*T2C.err*ts));%0.5<=qlye2<=10
%     T2C.err0 = T2C.err;
%     T2C.u0 = qlye(i,2);
%     
%     %PI controller qlye3-T3
%     qlye(i,3) = min(10e3,max(0.5e3,T3C.u0 + T3C.Kc*(T3C.err-T3C.err0) + T3C.Ki*T3C.err*ts));%0.5<=qlye3<=10
%     T3C.err0 = T3C.err;
%     T3C.u0 = qlye(i,3);
%     
% %             %PI controller for Vel-Telin
% %             VEl(i) = TElinC.u0 + TElinC.Kc*(TElinC.err-TElinC.err0) + TElinC.Ki*TElinC.err*ts;
% %             TElinC.err0 = TElinC.err;
% %             TElinC.u0 = VEl(i);
% %             V_El(i,:) = TElinC.u0*ones(1,par.N);
%     
%     %Proportional controller qcw-Telin
%     q_cw(i) = min(80e3,max(1e3,TElinC.u0 + TElinC.Kc*(TElinC.err)));
%     TElinC.u0 = q_cw(i);
     
    %% Region 2; 5.1<=P_inp<5.84
        T2C.err = T2C.set - Tk(i,2);
        T3C.err = T3C.set - Tk(i,3);
        PC.err = PC.set - Pcons(i);
        Iden1C.err = Iden1C.set - Iden1(i);
        TElinC.err = TElinC.set - T_ElinC(i);
    %     SOC.err = SOC.c_set - SOC.c(i);
    
        %PI controller qlye1-Iden1
        qlye(i,1) = min(10e3,max(0.5e3,Iden1C.u0 + Iden1C.Kc*(Iden1C.err-Iden1C.err0) + Iden1C.Ki*Iden1C.err*ts));%0.5<=qlye1<=10
        Iden1C.err0 = Iden1C.err;
        Iden1C.u0 = qlye(i,1);
    
        %PI controller qlye2-T2
        qlye(i,2) = min(10e3,max(0.5e3,T2C.u0 + T2C.Kc*(T2C.err-T2C.err0) + T2C.Ki*T2C.err*ts));%0.5<=qlye2<=10
        T2C.err0 = T2C.err;
        T2C.u0 = qlye(i,2);
    
        %PI controller qlye3-T3
        qlye(i,3) = min(10e3,max(0.5e3,T3C.u0 + T3C.Kc*(T3C.err-T3C.err0) + T3C.Ki*T3C.err*ts));%0.5<=qlye3<=10
        T3C.err0 = T3C.err;
        T3C.u0 = qlye(i,3);
    
        %Proportional controller Vel-P_inp
        VEl(i) = PC.u0 + PC.Kc*(PC.err);
        PC.u0 = VEl(i);
        V_El(i,:) = PC.u0*ones(1,par.N);
    
    
        %     %PI controller for qcw-SOC
        %     q_cw(i) = min(80e3,max(1e3,SOC.u0 + SOC.Kc*(SOC.err-SOC.err0) + SOC.Ki*SOC.err*ts));
        %     SOC.err0 = SOC.err;
        %     SOC.u0 = q_cw(i);

        %P controller for qcw-Telin
        q_cw(i) = min(80e3,max(1e3,TElinC.u0 + TElinC.Kc*(TElinC.err)));
        TElinC.u0 = q_cw(i);
        
    %% Region 3; 3.4<=P_inp<5.1
    
%         T1C.err = T1C.set - Tk(i,1);
%         T2C.err = T2C.set - Tk(i,2);
%         T3C.err = T3C.set - Tk(i,3);
%         PC.err = PC.set - Pcons(i);
%         TElinC.err = TElinC.set - T_ElinC(i);
%         %SOC.err = SOC.c_set - SOC.c(i);
%     
%             %PI controller qlye1-T1
%             qlye(i,1) = min(10e3,max(0.5e3,T1C.u0 + T1C.Kc*(T1C.err-T1C.err0) + T1C.Ki*T1C.err*ts));%0.5<=qlye1<=10
%             T1C.err0 = T1C.err;
%             T1C.u0 = qlye(i,1);
%     
%             %     T1C.u = q_lyek(1);
%             %     Iden1C.u = q_lyek(1);
%             %
%             %         %PI controller qlye1-T1/Iden1
%             %         T1C.u = min(10e3,max(0.5e3,T1C.u0 + T1C.Kc*(T1C.err-T1C.err0) + T1C.Ki*T1C.err*ts));%0.5<=qlye1<=10
%             %         T1C.err0 = T1C.err;
%             %         T1C.u0 = qlye(i,1);
%             %
%             %
%             %         Iden1C.u = min(10e3,max(0.5e3,Iden1C.u0 + Iden1C.Kc*(Iden1C.err-Iden1C.err0) + Iden1C.Ki*Iden1C.err*ts));%0.5<=qlye1<=10
%             %         Iden1C.err0 = Iden1C.err;
%             %         Iden1C.u0 = qlye(i,1);
%             %
%             %         if Iden1(i) < Iden1C.set
%             %             qlye(i,1) =  T1C.u;
%             %         else
%             %             qlye(i,1) =  max(Iden1C.u,T1C.u);
%             %         end
%     
%             %PI controller qlye2-T2
%             qlye(i,2) = min(10e3,max(0.5e3,T2C.u0 + T2C.Kc*(T2C.err-T2C.err0) + T2C.Ki*T2C.err*ts));%0.5<=qlye2<=10
%             T2C.err0 = T2C.err;
%             T2C.u0 = qlye(i,2);
%     
%             %PI controller qlye3-T3
%             qlye(i,3) = min(10e3,max(0.5e3,T3C.u0 + T3C.Kc*(T3C.err-T3C.err0) + T3C.Ki*T3C.err*ts));%0.5<=qlye3<=10
%             T3C.err0 = T3C.err;
%             T3C.u0 = qlye(i,3);
%     
%             %Proportional controller Vel-P_inp
%             VEl(i) = PC.u0 + PC.Kc*(PC.err);
%             PC.u0 = VEl(i);
%             V_El(i,:) = PC.u0*ones(1,par.N);
%     
%             %PI controller for qcw-Telin
% %             q_cw(i) = min(80e3,max(1e-2,TelinC.u0 + TelinC.Kc*(TelinC.err-TelinC.err0) + TelinC.Ki*TelinC.err*ts));
% %             TelinC.err0 = TelinC.err;
% %             TelinC.u0 = q_cw(i);
%             
%             %P controller for qcw-Telin 
%             q_cw(i) = min(80e3,max(1e3,TElinC.u0 + TElinC.Kc*(TElinC.err)));
%             TElinC.u0 = q_cw(i);
    
    %% Plant simulation
    %i = timestamp
    %j = electrolyzer sequence
    r =  F('x0',x0,'z0',z0,'p',[V_El(i,:), qlye(i,:), q_cw(i), pstoH2set(i), pstoO2set(i), Mass_Btset(i)]);
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
        U(i,j) = full(r.zf(j));                 %voltage/cell, [V]
        I(i,j) = full(r.zf(N+j));               %current, [A]
        P(i,j) = full(r.zf(2*N+j));             %power, [Watts]
        nH2elout(i,j) = full(r.zf(4*N+j));      %hydrogen production rate from individual electrolyzer, [mol/s]
        
        Temp(i,j)= full(r.xf(j));               %temperature of electrolyzers at all timestamps, [celsius]
        
        I_den(i,j) = 0.1*I(i,j)/par.EL(j).A;    %current density, [mA/cm^2]
        
        Qloss(i,j) = par.TherMo(j).A_El*(par.TherMo(j).hc*(Temp(i,j)-par.EL(j).Ta) + ...
            par.sigma*par.em*((Temp(i,j)+273.15)^4-(par.EL(j).Ta+273.15)^4));
        %heat loss to surrounding in the electrolyzer, [watts]
        
        Qgen(i,j) = par.EL(j).nc*(U(i,j)-par.EL(j).Utn)*I(i,j);
        %heat generated in the electrolyzer, [watts]
        
        Qlosslye(i,j) = qlye(i,j)*par.Const.CpLye*(Telin(i)-Temp(i,j));
        %heat taken out by the lye from the electrolyzer, [watts]
        
        Qnet(i,j) = Qgen(i,j)+Qlosslye(i,j)-Qloss(i,j);
        
        V_H2(i,j) = nH2elout(i,j)*0.0224136*3600;   %hydrogen production rate from individual electrolyzer, [Nm3/h]
        Ps(i,j) = P(i,j)/(1000*V_H2(i,j));          %Specific electricity consumption, [kWh/Nm3]
    end
    
    
    %% Calculate the input trajectory with PI controller
    % for storage pressure of the hydrogen tank
    Kc_pstoH2PI = Kc(1);                %controller gain
    taui_pstoH2PI = tauI(1);                %integral time constant
    e_pstoH2 = (PstoH2(i) - pstoH2set(i));%error term
    ZH2(i) = PIcontroller(zH2,Kc_pstoH2PI,taui_pstoH2PI,e_pstoH2,eint_pstoH2(i));
    
    
    %for storage pressure of the oxygen tank
    Kc_pstoO2PI = Kc(2);                %controller gain
    taui_pstoO2PI = tauI(2);                %integral time constant
    e_pstoO2 = (PstoO2(i) - pstoO2set(i));
    ZO2(i) = PIcontroller(zO2,Kc_pstoO2PI,taui_pstoO2PI,e_pstoO2,eint_pstoO2(i));
    
    
    %for mass in the buffer tank
    Kc_MassBtPI = Kc(3);
    taui_MassBtPI = tauI(3);
    e_Mbt = mBufferT(i) - Mass_Btset(i);
    q_H2O(i) = PIcontroller(Qwater,Kc_MassBtPI,taui_MassBtPI,e_Mbt,eint_Mbt(i));
    
    P_net(i)=sum(P(i,:));
    
    if rem(i,1000)==0
        disp(i)
    end
    
    Tk(i+1,:) = Temp(i,:);
    Pcons(i+1) = P_net(i);
    Iden1(i+1) = I_den(i,1);
    T_ElinC(i+1) = Telin(i);
% T_ElinC(i+1) = max(Temp(i,:))-Telin(i);
    %     SOC.ymeas = [Telin(i) Tw_out(i)]';
    %     SOC.c(i+1) = SOC.H*SOC.ymeas;
end

%% Plotting the results

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
% subplot(2,1,1)
% plot(q_cw)
% xlabel('Time, s')
% ylabel('q_{cw}, g/s')
% subplot(2,1,2)
% plot(SOC.c)
% xlabel('Time, s')
% ylabel('SOC')


figure()
subplot(321)
plot(qlye)
xlabel('Time, s')
ylabel('Lye flowrate El 2, g/s')
subplot(322)
plot(Temp)
xlabel('Time, s')
ylabel('T_k')
subplot(323)
plot(q_cw./1000)
xlabel('Time, s')
ylabel('q_{cw}, kg/s')
subplot(324)
plot(Telin)
xlabel('Time, s')
ylabel('T_{El_{in}}')
subplot(325)
plot(V_El)
xlabel('Time, s')
ylabel('Voltage')
subplot(326)
plot(I_den)
xlabel('Time, s')
ylabel('I_{den}, mA/cm^2')

%% Creating the data file
% save('data_closeloop_qlye1step_40hr5MW')
% load data_q1step3000_MVQcool_12hr
