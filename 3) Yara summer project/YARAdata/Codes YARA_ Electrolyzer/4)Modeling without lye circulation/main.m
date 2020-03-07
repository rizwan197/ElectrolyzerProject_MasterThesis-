clc
clear
close all

%% Load CasADi
addpath('/Users/mdrizwan/Documents/MATLAB/casadi-osx-matlabR2015a-v3.5.1')
import casadi.*


%Using CasADi we are solving system of ODE and nonlinear algebric eqns simultaneously
%1)UI*nc-Power = 0;
%2)U - (((r1+r2*T)*I)/A) - s*log10(((t1+(t2/T)+(t3/T^2))*I/A)+1) - Urev = 0;
%3)dT/dt = Qgen - Qloss - Qcool;
%4)Feff - ((.1*I/A)^2)/(f1+((.1*I/A)^2))*f2
%5)nH2 - Feff*nc*I/(ze*FC)
%6)nH2out - kvlvH2*VdispH2*sqrt(PstoH2-PoutH2)
%7)dPstoH2/dt = (TstoH2*Rg/VstoH2)*(nH2-nH2out)
%8)nO2out - kvlvO2*VdispO2*sqrt(PstoO2-PoutO2)
%9)dPstoO2/dt = (TstoO2*Rg/VstoO2)*(nO2-nO2out)

%% Loading parameters
global N
N = 3;                               %no. of electrolyzers
par = parElectrolyzer(N);

ze = par.Const.ze;                   %number of electrons transferred per reaction
FC = par.Const.FC;                   %faraday constant
R = par.Const.R;                     %gas constant, [J mol^-1 K^-1,]
alpha = par.Comp.alpha;              %compressor efficiency, 63%
k = par.Comp.k;                      %polytropic exponent, 1.62
Tel = par.Comp.Tel;                  %temperature of inlet to compressor after cooler, [Kelvin]
Pel = par.Comp.Pel;                  %pressure of electrolyzer outlet, [bar]
VstoH2 = par.Storage.VstoH2;         %volume of H2 storage, litres
VstoO2 = par.Storage.VstoO2;         %volume of O2 storage,litres
PoutH2 = par.Storage.PoutH2;         %outlet H2 pressure, bar
PoutO2 = par.Storage.PoutO2;         %outlet O2 pressure, bar
TstoH2 = par.Storage.TstoH2;         %temp of H2 storage, kelvin
TstoO2 = par.Storage.TstoO2;         %temp of H2 storage, kelvin
Rg = par.Storage.Rg;                 %gas constant in l bar K^-1 mol^-1
VdispH2 = par.Storage.VdispH2;       %steady state valve displacement for hydrogen outlet
VdispO2 = par.Storage.VdispO2;       %steady state valve displacement for oxygen outlet

%% Inputs for the simulation
Power = 21000;                          %net power to the plant, [Watts]
num_hr = .25;                             %no. of hours
t0 = 1;                                 %start, [s)]
ts = 1;                                 %time step, [s]
tf = num_hr*60*60;                      %final, [s]
tsamp = t0:ts:tf;
len = length(tsamp);                    %number of simulation time steps
%% Steady state solution and calculation of valve constant
%The steady state solution is calculated independently for each individual
%electrolyzer and then used as initial state for solving the dynamic model

%Valve constant is calculated for steady state condition and is kept
%constant for rest of the simulation

T_ini = 80; %assumed steady state temperature

nH2ss = zeros(N,1);     %H2 flow rate from electrolyzer at steady state, [mol/s]
nO2ss = zeros(N,1);     %O2 flow rate from electrolyzer at steady state, [mol/s]
nH2sout = zeros(N,1);   %H2 flow rate from storage at steady state, [mol/s]
nO2sout = zeros(N,1);   %O2 flow rate from storage at steady state, [mol/s]
Qloss = zeros(N,1);     %heat loss to ambient from electrolyzer, [Watts]
Qgen = zeros(N,1);      %internal heat generation, [Watts]

z00 = zeros(1,6*N);     %variable to store steady state solution of all electrolyzers
k_num=1;                %counter for creating steady state solution vector of all electrolyzers

param = struct([]);
for j=1:N %sequence from 1 to N electrolyzer
    
    %parameter assignmemt
    nc = par.EL(j).nc;    
    
    %initial guess for solving steady state
    param(j).Pel = Power;                       %power input to jth electrolyzer, [Watts]
    %    param(1).Pel = Power*.3;                   %option for assigning fraction of total power to jth electrolyser
    %    param(2).Pel = Power*.7;
    param(j).u0 = 1.3;                          %initial guess for voltage, [V]
    param(j).i0 = param(j).Pel/(nc*param(j).u0);%initial guess for current, [A]
    param(j).T0 = T_ini;
    
    %steady state solution
    [x0,z0] = elss(param,j);
    
    Qgen(j) = full(x0(1));      %rate of heat generation in jth electrolyzer at steady state
    Qloss(j) = full(x0(2));     %rate of heat loss in jth electrolyzer at steady state
    
    nH2ss(j) = full(z0(4));     %hydrogen flow rate from ith electrolyzer
    nO2ss(j) = full(z0(6));     %oxygen flow rate from ith electrolyzer
    nH2sout(j) = full(z0(5));   %contribution of jth electrolyzer to total hydrogen flowing out of the storage
    nO2sout(j) = full(z0(7));   %contribution of jth electrolyzer to total oxygen flowing out of the storage
    
    z00(k_num:j*6) = [full(z0(1:4)), full(x0)];%stacking of variables for initialization of dynamic eqns[U I Feff nH2in Qgen Qloss]
    k_num = k_num+6;
end

% Calculation of valve constants
nH2ss = sum(nH2ss);                 %net hydrogen flowrate from all electrolyzers at steady state(sum of individual contributions), [mol/s]
nO2ss = sum(nO2ss);                 %net oxygen flowrate from all electrolyzers at steady state(sum of individual contributions), [mol/s]
nH2sout = sum(nH2sout);             %net hydrogen flowrate from storage at steady state (sum of individual contributions), [mol/s]
nO2sout = sum(nO2sout);             %net oxygen flowrate from storage at steady state (sum of individual contributions), [mol/s]
kvalveH2 = nH2ss/VdispH2;           %valve constant for hydrogen outlet
kvalveO2 = nO2ss/VdispO2;           %valve constant for oxygen outlet

Psto_iniH2 = (nH2ss/(kvalveH2*VdispH2))^2 + PoutH2;         %initial H2 storage pressure (calculated from steady state solution) [bar]/can be set manually
Psto_iniO2 = (nO2ss/(kvalveO2*VdispO2))^2 + PoutO2;         %initial O2 storage pressure (calculated from steady state solution) [bar]/can be set manually

x0 = [T_ini Psto_iniH2 Psto_iniO2];                         %initial vector for dynamic solution of differential variables
z0 = [z00 nH2ss nH2sout nO2ss nO2sout];                     %initial vector for dynamic solution of algebriac variables

%calculation of Qcool at steady state
Qgen = sum(Qgen);                   %rate of heat generation in the electrolyzer stack at steady state
Qloss = sum(Qloss);                 %rate of heat loss in the electrolyzer stack at steady state
Qcool = Qgen-Qloss;                 %coolant flowrate at steady state

%% Manipulated variables
%these are the degree of freedoms that we will utilise to control the
%system

P = zeros(len,N);                            %power, [Watt], len is the length of time vector
for j=1:N
P(1:500,j)=param(j).Pel*1;               %incremental step change in power for all electrolysers
P(501:end,j)=param(j).Pel*1;
%P(1001:end,j)=1.2*param(j).Pel;
%P(1501:end,j)=0.8*param(j).Pel;
%P(2001:end,j)=param(j).Pel*1;
end
% P(301:end,2)=param(2).Pel*1.1;              %incremental step change in power for individual electrolysers
% P(1001:end,3)=param(3).Pel*0.95;

Qc = Qcool*ones(len,1);             %coolant flow rate as a manipulated variable
%Qc(1:500)=Qcool*1.2;               %incremental step change in cooling rate
%Qc(501:end)=Qcool*1.5;
%Qc(1201:end)=1.025*Qcool;

ZH2 = VdispH2*ones(len,1);          %H2 valve displacement as a manipulated variable
%ZH2(1201:end) = .1;                 %change in H2 valve displacement

ZO2 = VdispO2*ones(len,1);          %O2 valve displacement as a manipulated variable
%ZO2(501:end) = .7;                 %change in O2 valve displacement

%% Initialize plotting variables
Temp = zeros(len+1,1);                  %temp of the electrolyzer, [C]
Temp(1) = full(T_ini);

PstoH2 = zeros(len+1,1);                %H2 storage pressure, [bar]
PstoH2(1) = full(Psto_iniH2);

PstoO2 = zeros(len+1,1);                %O2 storage pressure, [bar]
PstoO2(1) = full(Psto_iniO2);

U = zeros(len,1);                       %voltage/cell in each of the electrolyzer, [V] 
I = zeros(len,1);                       %current in each electrolyzer, [A]
I_den = zeros(len,1);                   %current density in the electrolyzer, [A/m^4]
nH2in = zeros(len,1);                   %net hydrogen flow rate in to the storage, [mol/s]
nH2out = zeros(len,1);                  %net hydrogen flowrate out from the storage, [mol/s]
nH2elout = zeros(len,1);                %hydrogen flowrate from each of the individual electrolyzer, [mol/s]
nO2in = zeros(len,1);                   %net oxygen flow rate in to the storage, [mol/s]
nO2out = zeros(len,1);                  %net oxygen flowrate out from the storage, [mol/s]
PcompH2 = zeros(len,1);                 %compressor power for hydrogen, [watts]
PcompO2 = zeros(len,1);                 %compressor power for oxygen, [watts]

%% Solving for T_el,U,I,Feff,nH2,nO2

%T=x(1);PstoH2=x(2);PstoO2=x(3);
%U=z(6j-5);I=z(6j-4);Feff=z(6j-3);nH2=z(6j-2);Qgen=z(6j-1);Qloss=z(6j);nH2out=z(6N+2);nO2=z(6N+3);nO2out=z(6N+4)

eqnAlg = SX.zeros(6*N+4,1);
eqnDiff = SX.zeros(3,1);
z = SX.sym('z',6*N+4); x = SX.sym('x',3); p = SX.sym('p',N+3);
%standard casadi notation, z: algebric variable, x: differential variable,
%p: parameters (MV)

%here, we are now writing equations for transient behaviour of the electrolyzer. 
for j=1:N

%parameter assignment
r1 = par.U(j).r1;
r2 = par.U(j).r2;
s = par.U(j).s;
t1 = par.U(j).t1;
t2 = par.U(j).t2;
t3 = par.U(j).t3;
f1 = par.U(j).f1;
f2 = par.U(j).f2;

Ct = par.TherMo(j).Ct;
Rt = par.TherMo(j).Rt;

Utn = par.EL(j).Utn;
Urev = par.EL(j).Urev;
nc = par.EL(j).nc;
A = par.EL(j).A;
Ta = par.EL(j).Ta;
Tstd = par.EL(j).Tstd;

%system of algebric equations for an electrolyzer 
eqnAlg(6*j-5) = z(6*j-5)*z(6*j-4)*nc- p(j);                     %power = nc*UI
eqnAlg(6*j-4) = z(6*j-5) - (r1+r2*x(1))*z(6*j-4)/A - s*log10(((t1+t2/x(1)+t3/x(1)^2)*z(6*j-4)/A)+1) - Urev; %U-I relationship
eqnAlg(6*j-3) = z(6*j-3) - ((.1*z(6*j-4)/A)^2)/(f1+((.1*z(6*j-4)/A)^2))*f2;                                 %faraday efficiency
eqnAlg(6*j-2) = z(6*j-2) - z(6*j-3)*nc*z(6*j-4)/(ze*FC);        %nH2 from individual electrolyzer
eqnAlg(6*j-1) = z(6*j-1) - nc*(z(6*j-5)-Utn)*z(6*j-4);          %Qgen eqn for individual electrolyzer
eqnAlg(6*j) = z(6*j) - (x(1)-Ta)/Rt;                            %Qloss for individual electrolyzer


end
sum_H2net = SX.zeros(1,1);
Qgen_net = SX.zeros(1,1);
Qloss_net = SX.zeros(1,1);

for j=1:N
sum_H2net = sum_H2net + z(6*j-2);   %sum of hydrogen from all individual electrolyzers
Qgen_net = Qgen_net + z(6*j-1);     %total heat generated from all elctrolyzers 
Qloss_net = Qloss_net + z(6*j);     %total heat loss from all electrolyzers 
end

eqnAlg(6*N+1) = z(6*N+1) - sum_H2net;                           %algebric eqn for net hydrogen flowrate from all the electrolyzers
eqnAlg(6*N+2) = z(6*N+2) - kvalveH2*p(N+2)*sqrt(x(2)-PoutH2);   %algebric eqn for net hydrogen flowrate from the storage tank
eqnAlg(6*N+3) = z(6*N+3) - z(6*N+1)/2;                          %algebric eqn for net oxygen flowrate from all the electrolyzers
eqnAlg(6*N+4) = z(6*N+4) - kvalveO2*p(N+3)*sqrt(x(3)-PoutO2);   %algebric eqn for net oxygen flowrate from the storage tank

eqnDiff(1) = (Qgen_net - Qloss_net - p(N+1))/(Ct*N);            %differential eqn for the lye temperature
eqnDiff(2) = (TstoH2*Rg/VstoH2)*(z(6*N+1)-z(6*N+2));            %differential eqn for hydrogen storage pressure
eqnDiff(3) = (TstoO2*Rg/VstoO2)*(z(6*N+3)-z(6*N+4));            %differential eqn for oxygen storage pressure

dae = struct('x',x,'z',z,'p',p,'ode',eqnDiff,'alg',eqnAlg);
F = integrator('F', 'idas', dae);


for i=1:len
%i = timestamp
%j = electrolyzer sequence

r = F('x0',x0,'z0',z0,'p',[P(i,:) Qc(i) ZH2(i) ZO2(i)]);
x0 = full(r.xf);            %updating solution as new initial conditions
z0 = full(r.zf);

%% Storing values in plotting variables

%calculation of compressor power
PcompH2(i) = full(((r.zf(6*N+1)*k*R*Tel)/(alpha*(k-1)))*(((r.xf(2)/Pel)^((k-1)/k))-1));
PcompO2(i) = full(((r.zf(6*N+3)*k*R*Tel)/(alpha*(k-1)))*(((r.xf(3)/Pel)^((k-1)/k))-1));%assuming same k and Tel for O2



Temp(i+1) = full(r.xf(1));              %temperature at all timestamps, [celsius]
PstoH2(i+1) = full(r.xf(2));            %hydrogen storage pressure at all timestamps, [bar]
PstoO2(i+1) = full(r.xf(3));            %oxygen storage pressure at all timestamps, [bar]
nH2in(i) = full(r.zf(6*N+1));           %net hydrogen flow rate in to the storage at all timestamps, [mol/s]
nH2out(i) = full(r.zf(6*N+2));          %net hydrogen flowrate out from the storage at all timestamps, [mol/s]


nO2in(i) = full(r.zf(6*N+3));           %net oxygen flow rate in to the storage at all timestamps, [mol/s]
nO2out(i) = full(r.zf(6*N+4));          %net oxygen flowrate out from the storage at all timestamps, [mol/s]

    for j=1:N
        U(i,j) = full(r.zf(6*j-5));         %voltage/cell, [V]
        I(i,j) = full(r.zf(6*j-4));         %current, [A]
        I_den(i,j) = 0.1*I(i,j)/A;          %current density, [A/m^4]
        nH2elout(i,j) = full(r.zf(6*j-2));  %hydrogen flowrate from individual electrolyzer, [mol/s]
        nH2outflow(i,j) = (nH2out(i)*0.0224136*3600)/par.EL(j).nc; %net hydrogen flowrate out from the storage at all timestamps, [Nm3/h]
    end
    
    disp(i)
end

%% Plotting the results
figure()
subplot(2,1,2)
plot(I)
xlabel('Time, s')
ylabel('Current, A')
legend(strcat('n_c = ',num2str([15 21 27]')))
% subplot(2,1,2)
% plot(U)
% xlabel('Time, s')
% ylabel('Voltage, V')
% legend(strcat('Electrolyzer ',num2str([1 2 3]')))
subplot(2,1,1)
plot(nH2elout)
xlabel('Time, s')
ylabel('H_2 Flowrate, mol/s')
legend(strcat('El : ',num2str([1 2 3]')))
% subplot(3,1,3)
% plot(I_den,U,'b*')
% xlabel('Current Density, mA/cm^2')
% ylabel('V/cell')

figure()
% subplot(2,1,1)
plot(Temp)
xlabel('Time, s')
ylabel('Temperature, Celsius')
ylim([78 82])
% grid on
% subplot(2,1,2)
% plot(nH2elout)
% xlabel('Time, s')
% ylabel('H_2 Flowrate, mol/s')
% legend(strcat('Electrolyzer ',num2str([1 2 3]')))
% grid on

% plot(PstoH2)
% xlabel('Time, s')
% ylabel('H_2 Storage pressure, bar')
% subplot(3,1,3)
% plot(PstoO2)
% xlabel('Time, s')
% ylabel('O_2 Storage pressure, bar')

figure()
plot(nH2in)
xlabel('Time, s')
ylabel('Flowrate, mol/s')
hold on
plot (nH2out)
legend('n_H_2 in','n_H_2 out')
xlabel('Time, s')
hold off

figure()
plot(nH2outflow)
xlabel('Time, s')
ylabel('Flowrate, Nm3/h')

% figure()
% plot(nO2in)
% xlabel('Time, s')
% ylabel('Flowrate, mol/s')
% hold on
% plot (nO2out)
% legend('n_O_2 in','n_O_2 out')
% xlabel('Time, s')
% hold off

figure()
subplot(2,1,1)
plot(PstoH2)
xlabel('Time, s')
ylabel('H_2 Storage pressure, bar')
ylim([15,40])
grid on
subplot(2,1,2)
plot(PcompH2)
xlabel('Time, s')
ylabel('H_2 compressor Power, W')
grid on
% plot(PcompO2)
% xlabel('Time, s')
% ylabel('O_2 compressor Power, W')