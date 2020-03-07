clc
clear
close all

%% Load CasADi
addpath('/Users/mdrizwan/Documents/MATLAB/casadi-osx-matlabR2015a-v3.5.1')
import casadi.*

%Using CasADi we are solving system of ODE and nonlinear algebraic eqns simultaneously

%Nonlinear algebraic equation are:
%1N)UI*nc-Power = 0;
%2N)U - (((r1+r2*T)*I)/A) - s*log10(((t1+(t2/T)+(t3/T^2))*I/A)+1) - Urev = 0;
%3N)U*nc - V = 0; U=cell voltage; V=electrolyzer voltage
%4N)Feff - ((.1*I/A)^2)/(f1+((.1*I/A)^2))*f2;
%5N)nH2el - Feff*nc*I/(ze*FC);
%6N)qH2Oloss - nH2*MwtH2O = 0;
%7)nH2elnet - sum(nH2el) = 0;
%8)nH2out - kvlvH2*VdispH2*sqrt(PstoH2-PoutH2) = 0;
%9)nO2elnet - nH2elnet/2 = 0;
%10)nO2out - kvlvO2*VdispO2*sqrt(PstoO2-PoutO2) = 0;
%11)Tout - ((sum(qlye*T)*CpLye - sum(qloss*T)*Cp + sum(qloss)*(Cp-CpLye)*Tref)/((sum(qlye)-sum(qloss))*Cp) = 0;
%12)Tin - Tout + Qcool/qlye*CpLye = 0;

%ODE eqautions are:
%1N)dT/dt = qlye*CpLye*(T_in-T) + nc*(U-Utn)*I - (T-Ta)/Rt))/Ct;
%2)dPstoH2/dt = (TstoH2*Rg/VstoH2)*(nH2-nH2out);
%3)dPstoO2/dt = (TstoO2*Rg/VstoO2)*(nO2-nO2out);
%4)dM_bt/dt = (qlye-qloss) + qH2O - qlye;

%Parameters for the simulation are:
%1)Power
%2)qlye
%3)Qcool
%4)zH2
%5)zO2
%6)qH2O

%% Loading parameters
N = 3;                               %no. of electrolyzers
par = parElectrolyzer(N);
tau1 = 50;
%% Inputs for the simulation
Power = 2134000*3;                        %power to all the electrolyzers, [Watts]
num_hr = 1;                             %no. of hours
t0 = 1;                                 %start, [s)]
ts = 1;                                 %time step, [s]
tf = num_hr*60*60;                      %final, [s]
tsamp = t0:ts:tf;
len = length(tsamp);                    %number of simulation time steps
tstep = 200;
%% Steady state solution and calculation of valve constant
%The steady state solution is calculated by solving non linear algebraic equations for the steady state system
%Valve constant (k_vlv) is calculated for steady state condition and is kept constant for rest of the simulation

%initialization for steady state solution

%T_ini is regulatory layer CV i.e. CV2 but it is supervisory layer MV i.e.
%MV1
T_ini = 65*ones(len,1);                     %assumed steady state temperature of the lye into the electrolyzer
%T_ini(tstep:end,1) = 66;

param.P0 = Power/N;            %initial guess for power input to each electrolyzer, [Watts]
param.u0 = 1.5;                   %initial guess for cell voltage, [V]
param.T0 = T_ini(1);

q_lye = N*6648;              %total lye circulated, [g/s]
for nEl = 1:N
    param.q_lye(nEl) = q_lye/N; %lye through each electrolyzer, [g/s]
end

% z0s = [u_ini, i_ini, Pk, Feff_ini, nH2el_ini, qH2O_loss, nH2out_ini, nO2el_ini, nO2out_ini, T_out Tw_out]; (7N+4)
% x0s = [Tk_ini, Qcool, V_ini, q_cw, UA_hex];(2N+3)

[x0s,z0s] = elss(param,par);

z00=[];
for j=1:N
    z00 = [z00 z0s(j) z0s(N+j) z0s(2*N+j) z0s(3*N+j) z0s(4*N+j) z0s(5*N+j)];    %stacking of variables for initialization of dynamic eqns[U I P Feff nH2in qH2Oloss]
end

nH2ss = z0s(4*N+1:5*N);      %hydrogen flow rate from kth electrolyzer
nO2ss = z0s(6*N+2:7*N+1);    %oxygen flow rate from kth electrolyzer
nH2sout = z0s(6*N+1);        %total hydrogen flowing out of the storage at steady state
nO2sout = z0s(7*N+2);        %total oxygen flowing out of the storage at steady state

% Calculation of valve constants
nH2ss = sum(nH2ss);                             %net hydrogen flowrate from all electrolyzers at steady state(sum of individual contributions), [mol/s]
nO2ss = sum(nO2ss);                             %net oxygen flowrate from all electrolyzers at steady state(sum of individual contributions), [mol/s]
kvalveH2 = nH2ss/par.Storage.VdispH2;           %valve constant for hydrogen outlet
kvalveO2 = nO2ss/par.Storage.VdispO2;           %valve constant for oxygen outlet

Psto_iniH2 = (nH2ss/(kvalveH2*par.Storage.VdispH2))^2 + par.Storage.PoutH2; %initial H2 storage pressure (calculated from steady state solution) [bar]
Psto_iniO2 = (nO2ss/(kvalveO2*par.Storage.VdispO2))^2 + par.Storage.PoutO2; %initial O2 storage pressure (calculated from steady state solution) [bar]
mass_bt0 = 1000000;                                                            %mass of the liquid in the buffer tank at steady state, [g]

%initial values when using dynamic model for the heat exchanger 
% x0 = [x0s(1:N) Psto_iniH2 Psto_iniO2 mass_bt0 T_ini(1) z0s(7*N+4)];                             %initial vector for dynamic solution of differential variables
% z0 = [z00 nH2ss nH2sout nO2ss nO2sout z0s(7*N+3)];         %initial vector for dynamic solution of algebriac variables

%initial values when using steady state model for the heat exchanger
x0 = [x0s(1:N) Psto_iniH2 Psto_iniO2 mass_bt0 T_ini(1)];                             %initial vector for dynamic solution of differential variables
z0 = [z00 nH2ss nH2sout nO2ss nO2sout z0s(7*N+3) T_ini(1) z0s(7*N+4)];         %initial vector for dynamic solution of algebriac variables

%defining initial value of parameters for the simualation
qf_cw = x0s(2*N+2);                     %cooling water flow rate, [g/s]
UA_hex = x0s(2*N+3);                    %UA of the heat exchanger, [Watt/K]
Vss = x0s(N+2:2*N+1);                   %voltage across electrolyzers, [Volts]
Qwater = sum(z0s(5*N+1:6*N));           %total water lost during electrolysis, [grams/sec]

%% Manipulated variables
%these are the degree of freedoms that we will utilise to control the system

V_El = zeros(len,N);                       %voltage across the electrolyzer, [Watt], len is the length of time vector
for j = 1:N
    V_El(1:end,j) = Vss(j)*1;               %incremental step change in common voltage across all electrolysers
    V_El(tstep:end,j)=Vss(j)*1;
end

qlye = zeros(len,N);                        %lye flowrate, [g/s]
for j = 1:N
    qlye(1:end,j) = param.q_lye(j)*1;       %assumed same lye flowarate to all the electrolyzers
end
qlye(tstep:end,1) = param.q_lye(j)*1;

q_cw = qf_cw*ones(len,1);                     %cooler duty as a manipulated variable, [J/s]
q_cw(tstep:end) = qf_cw*1.2;                  %incremental step change in cooling water flowrate

ZH2 = par.Storage.VdispH2*ones(len,1);      %H2 valve displacement as a manipulated variable
%ZH2(tstep:end) = .2;                         %change in H2 valve displacement

ZO2 = par.Storage.VdispO2*ones(len,1);      %O2 valve displacement as a manipulated variable
%ZO2(tstep:end) = .7;                         %change in O2 valve displacement

qH2O = Qwater*ones(len,1);                  %flow rate of water added to buffer tank as a manipulated variable, [g/s]
%qH2O(tstep:end)=Qwater*1.2;                  %incremental step change in the water flow rate

%% Initialize plotting variables
Temp = zeros(len+1,N);                  %temp of the electrolyzer, [C]
Temp(1,:) = full(T_ini(1));

PstoH2 = zeros(len+1,1);                %H2 storage pressure, [bar]
PstoH2(1) = full(Psto_iniH2);

PstoO2 = zeros(len+1,1);                %O2 storage pressure, [bar]
PstoO2(1) = full(Psto_iniO2);

U = zeros(len,1);                       %voltage/cell in each of the electrolyzer, [V]
I = zeros(len,1);                       %current in each electrolyzer, [A]
P = zeros(len,1);
I_den = zeros(len,1);                   %current density in the electrolyzer, [A/m^2]
nH2in = zeros(len,1);                   %net hydrogen flow rate in to the storage, [mol/s]
nH2out = zeros(len,1);                  %net hydrogen flowrate out from the storage, [mol/s]
nH2elout = zeros(len,1);                %hydrogen flowrate from each of the individual electrolyzer, [mol/s]
nO2in = zeros(len,1);                   %net oxygen flow rate in to the storage, [mol/s]
nO2out = zeros(len,1);                  %net oxygen flowrate out from the storage, [mol/s]
Tout = zeros(len,1);
Tin = zeros(len,1);
level = zeros(len,1);
SpecEl = zeros(len,1);
PcompH2 = zeros(len,1);                 %compressor power for hydrogen, [watts]
PcompO2 = zeros(len,1);                 %compressor power for oxygen, [watts]
Qloss = zeros(len,1);                   %heat loss to surrounding in the electrolyzer, [watts]
Qgen = zeros(len,1);                    %heat generated in the electrolyzer, [watts]
Qlosslye = zeros(len,1);                %heat taken out by the lye from the electrolyzer, [watts]
nH2inSto = zeros(len,1);                %net hydrogen flowrate into the storage at all timestamps, [Nm3/h]
P_net = zeros(len,1);                   %net power to the electrolyzer assembly, [watts]

%% Solving for T_el,U,I,Feff,nH2,nO2,Tout,Tin and P_sto and V_bt

%differential variables:
%T=x(1:N);PstoH2=x(N+1);PstoO2=x(N+2);V_bt=x(N+3);Tin=x(N+4);Tw_out=x(N+5)
%algebraic varibales: %U=z(6j-5);I=z(6j-4);P=z(6j-3);Feff=z(6j-2);nH2=z(6j-1);qH2Oloss=z(6j);nH2net=z(6N+1);nH2out=z(6N+2);nO2net=z(6N+3);nO2out=z(6N+4);
%Tout=z(6N+5);
%parameters: V_El=p(1:N);qlye_k=p(N+1:2N);q_cw=p(2*N+1);zH2=p(2*N+2);zO2=p(2*N+3);qH2O=p(2*N+4)

eqnAlg = SX.zeros(6*N+7,1);
eqnDiff = SX.zeros(N+4,1);
z = SX.sym('z',6*N+7); x = SX.sym('x',N+4); p = SX.sym('p',2*N+4);
%standard casadi notation, z: algebraic variable, x: differential variable,
%p: parameters (MV)

for j = 1:N
    
    %T=x(1:N);PstoH2=x(N+1);PstoO2=x(N+2);V_bt=x(N+3);Tin=x(N+4);Tw_out=x(N+5)
    %algebraic varibales: %U=z(6j-5);I=z(6j-4);P=z(6j-3);Feff=z(6j-2);nH2=z(6j-1);qH2Oloss=z(6j);nH2net=z(6N+1);nH2out=z(6N+2);nO2net=z(6N+3);nO2out=z(6N+4);
    %Tout=z(6N+5);
    %parameters: V_El=p(1:N);qlye_k=p(N+1:2N);q_cw=p(2*N+1);zH2=p(2*N+2);zO2=p(2*N+3);qH2O=p(2*N+4)
    
    %calculation of Urev as a function of temperature
    par.EL(j).Urev = 1.5184 - 1.5421e-3*(273+x(j)) + 9.523e-5*(273+x(j))*log((273+x(j))) + ...
        9.84e-8*(273+x(j))^2; %relation for Urev with T from LeRoy eqn. 58
    
    %system of algebraic equations for an electrolyzer
    eqnAlg(6*j-5) = z(6*j-5)*z(6*j-4)*par.EL(j).nc - z(6*j-3);                                                            %power = nc*UI
    eqnAlg(6*j-4) = z(6*j-5) - (par.U(j).r1+par.U(j).r2*x(j))*z(6*j-4)/par.EL(j).A - par.U(j).s*log10(((par.U(j).t1+par.U(j).t2/x(j)+...
        par.U(j).t3/x(j)^2)*z(6*j-4)/par.EL(j).A)+1) - par.EL(j).Urev;                                                    %U-I relationship
    eqnAlg(6*j-3) = z(6*j-5)*par.EL(j).nc - p(j);                                                                         %U(j).nc(j)=V
    eqnAlg(6*j-2) = z(6*j-2) - ((.1*z(6*j-4)/par.EL(j).A)^2)/(par.U(j).f1+((.1*z(6*j-4)/par.EL(j).A)^2))*par.U(j).f2;     %faraday efficiency
    eqnAlg(6*j-1) = z(6*j-1) - z(6*j-2)*par.EL(j).nc*z(6*j-4)/(par.Const.ze*par.Const.FC);                                %nH2, H2 production rate from individual electrolyzer
    eqnAlg(6*j) = z(6*j) - z(6*j-1)*par.Const.Mwt;                                                                        %flowrate of water lost, [g/s]
end

sum_H2net = SX.zeros(1,1);
netqlye = SX.zeros(1,1);
netqloss = SX.zeros(1,1);
qlyeT = SX.zeros(1,1);
qlossT = SX.zeros(1,1);

for j = 1:N
    sum_H2net = sum_H2net + z(6*j-1);       %sum of hydrogen from all individual electrolyzers
    netqlye = netqlye + p(N+j);             %sum the lye flowing into all individual electrolyzers
    netqloss = netqloss + z(6*j);           %sum of water lost from all individual electrolyzers
    qlyeT = qlyeT+p(N+j)*x(j);              %calculate term qlye(k).T(k)
    qlossT = qlossT+z(6*j)*x(j);            %calculate term qH2Oloss(k).T(k)
end

% %epsilon-NTU method for lye circulation loop heat exchanger model
% Ntu = UA_hex/(p(2*N+1)*par.Const.Cp);
% C = (p(2*N+1)*par.Const.Cp)/(netqlye*par.Const.CpLye);
% epsC = (1-exp(-Ntu*(C-1)))/(C-exp(-Ntu*(C-1)));
% epsH = epsC*C;
% eqnAlg(6*N+6) = z(6*N+6) - ((1-epsH)*z(6*N+5) + epsH*par.Tw_in);                %calculation of T_in using epsilon-NTU approach, exit temperature of hot stream 
% eqnAlg(6*N+7) = z(6*N+7) - (epsC*z(6*N+5) + (1-epsC)*par.Tw_in);                %calculation of Tw_out using epsilon-NTU approach, exit temperature of the cold stream 

%LMTD approach
deltaT_LMTD = ((z(6*N+6)-par.Tw_in)-(z(6*N+5)-z(6*N+7)))/log((z(6*N+6)-par.Tw_in)/(z(6*N+5)-z(6*N+7)));
eqnAlg(6*N+6) = netqlye*par.Const.CpLye*(z(6*N+5)-z(6*N+6)) - (UA_hex*deltaT_LMTD);  %calculation of T_in using LMTD approach, exit temperature of hot stream 
eqnAlg(6*N+7) = p(2*N+1)*par.Const.Cp*(z(6*N+7)-par.Tw_in) - (UA_hex*deltaT_LMTD);   %calculation of Tw_out using LMTD approach, exit temperature of the cold stream

eqnAlg(6*N+1) = z(6*N+1) - sum_H2net;                                           %algebraic eqn for net hydrogen flowrate from all the electrolyzers
eqnAlg(6*N+2) = z(6*N+2) - kvalveH2*p(2*N+2)*sqrt(x(N+1)-par.Storage.PoutH2);   %algebraic eqn for net hydrogen flowrate from the storage tank
eqnAlg(6*N+3) = z(6*N+3) - z(6*N+1)/2;                                          %algebraic eqn for net oxygen flowrate from all the electrolyzers
eqnAlg(6*N+4) = z(6*N+4) - kvalveO2*p(2*N+3)*sqrt(x(N+2)-par.Storage.PoutO2);   %algebraic eqn for net oxygen flowrate from the storage tank
eqnAlg(6*N+5) = z(6*N+5) - (((qlyeT*par.Const.CpLye) - (qlossT*par.Const.Cp) + netqloss*(par.Const.Cp-par.Const.CpLye)*par.Const.Tref)/...
    ((netqlye-netqloss)*par.Const.CpLye));                                      %calculation of Tout


for j = 1:N
    eqnDiff(j) = (p(N+j)*par.Const.CpLye*(x(N+4)-x(j)) + par.EL(j).nc*(z(6*j-5)-par.EL(j).Utn)*z(6*j-4) - ...
        par.TherMo(j).A_surf*z(6*j-4)* par.EL(j).nc/1000*(par.TherMo(j).hc*(x(j)-par.EL(j).Ta) + ...
        par.sigma*par.em*((x(j)+273.15)^4-(par.EL(j).Ta+273.15)^4)))/(par.TherMo(j).CtS*z(6*j-4));%differential eqn for the electrolyzer temperature
end

eqnDiff(N+1) = (par.Storage.TstoH2*par.Storage.Rg/par.Storage.VstoH2)*(z(6*N+1)-z(6*N+2));  %differential eqn for hydrogen storage pressure
eqnDiff(N+2) = (par.Storage.TstoO2*par.Storage.Rg/par.Storage.VstoO2)*(z(6*N+3)-z(6*N+4));  %differential eqn for oxygen storage pressure
eqnDiff(N+3) = (netqlye-netqloss) + p(2*N+4) - netqlye;                                     %differential eqn for mass in the buffer tank, [grams i.e. pho*V]
eqnDiff(N+4) = (z(6*N+6) - x(N+4))/tau1;                                                        %first order filter to the Th_out to account for the dynamics

%dynamic thermal balance for the heat exchanger, T_in and Tw_out are
%differential variables now
% deltaT1 = x(N+4) - par.Tw_in;%difference in temperatures of hot and cold streams at inlet
% deltaT2 = z(6*N+5) - x(N+5);%difference in temperatures of hot and cold streams at outlet
% deltaT_LMTD = ((x(N+4)-par.Tw_in)-(z(6*N+5)-x(N+5)))/log((x(N+4)-par.Tw_in)/(z(6*N+5)-x(N+5)));
% deltaT_LMTD = (((deltaT1^(1/3))+(deltaT2)^(1/3))/2)^1/3;%deltaT_LMTD using underwood expression
% deltaT_LMTD = ((((deltaT1)+(deltaT2))/2)^(1/3))*((sqrt(deltaT1*deltaT2))^(2/3));%deltaT_LMTD 
% eqnDiff(N+4) = (z(6*N+5)-x(N+4))/tau1 - (UA_hex*deltaT_LMTD/(par.Const.rho*par.Const.CpLye*par.Const.Vh));%differential eqn for the hot stream exit temp from heat exchanger
% eqnDiff(N+5) = (par.Tw_in-x(N+5))/tau1 - (UA_hex*deltaT_LMTD/(par.Const.rho*par.Const.Cp*par.Const.Vc));%differential eqn for the cold stream exit temperature from heat exchanger

dae = struct('x',x,'z',z,'p',p,'ode',eqnDiff,'alg',eqnAlg);
F = integrator('F', 'idas', dae);


for i=1:len
    %i = timestamp
    %j = electrolyzer sequence
    
    r = F('x0',x0,'z0',z0,'p',[V_El(i,:) qlye(i,:) q_cw(i) ZH2(i) ZO2(i) qH2O(i)]);
    x0 = full(r.xf);            %updating solution as new initial conditions
    z0 = full(r.zf);
    
    
    %% Storing values in plotting variables
    %T=x(1:N);PstoH2=x(N+1);PstoO2=x(N+2);V_bt=x(N+3);Tin=x(N+4);Tw_out=x(N+5)
    %algebraic varibales: %U=z(6j-5);I=z(6j-4);P=z(6j-3);Feff=z(6j-2);nH2=z(6j-1);qH2Oloss=z(6j);
    %nH2net=z(6N+1);nH2out=z(6N+2);nO2net=z(6N+3);nO2out=z(6N+4);Tout=z(6N+5);
    %parameters: V_El=p(1:N);qlye_k=p(N+1:2N);q_cw=p(2*N+1);zH2=p(2*N+2);zO2=p(2*N+3);qH2O=p(2*N+4)
    
    %calculation of compressor power
    PcompH2(i) = full(((r.zf(6*N+1)*par.Comp.k*par.Const.R*par.Comp.Tel)/(par.Comp.alpha*(par.Comp.k-1)))*(((r.xf(N+1)/par.Comp.Pel)^((par.Comp.k-1)/par.Comp.k))-1));
    PcompO2(i) = full(((r.zf(6*N+3)*par.Comp.k*par.Const.R*par.Comp.Tel)/(par.Comp.alpha*(par.Comp.k-1)))*(((r.xf(N+2)/par.Comp.Pel)^((par.Comp.k-1)/par.Comp.k))-1));
    %assuming same k and Tel for O2
    
    PstoH2(i+1) = full(r.xf(N+1));          %hydrogen storage pressure at all timestamps, [bar]
    PstoO2(i+1) = full(r.xf(N+2));          %oxygen storage pressure at all timestamps, [bar]
    nH2in(i) = full(r.zf(6*N+1));           %net hydrogen flow rate in to the storage at all timestamps, [mol/s]
    nH2out(i) = full(r.zf(6*N+2));          %net hydrogen flowrate out from the storage at all timestamps, [mol/s]
    nO2in(i) = full(r.zf(6*N+3));           %net oxygen flow rate in to the storage at all timestamps, [mol/s]
    nO2out(i) = full(r.zf(6*N+4));          %net oxygen flowrate out from the storage at all timestamps, [mol/s]
    Tout(i) = full(r.zf(6*N+5));            %temperature of lye after mixing before going to the buffer tank, [celsius]
    Tin(i) = full(r.xf(N+4));             %temperature of lye going into the electrolyzer, [celsius]
    Tw_out(i) = full(r.zf(6*N+7));          %exit temperature of the cooling water, [celsius]
    level(i) = full(r.xf(N+3));
    
    for j=1:N
        U(i,j) = full(r.zf(6*j-5));             %voltage/cell, [V]
        I(i,j) = full(r.zf(6*j-4));             %current, [A]
        P(i,j) = full(r.zf(6*j-3));             %power, [Watts]
        Temp(i+1,j)= full(r.xf(j));             %temperature of electrolyzers at all timestamps, [celsius]
        I_den(i,j) = 0.1*I(i,j)/par.EL(j).A;    %current density, [mA/cm^2]
        nH2elout(i,j) = full(r.zf(6*j-1));      %hydrogen production rate from individual electrolyzer, [mol/s]
        SpecEl(i,j) = (P(i,j)*(10^-6))./(nH2elout(i,j)*par.Const.MwtH2*(10^-6)*3600);   %specific electricity, [MWh/tonne H2]
        Qloss(i,j) = par.TherMo(j).A_surf*I(i,j)* par.EL(j).nc/1000*(par.TherMo(j).hc*(Temp(i+1,j)-par.EL(j).Ta) + ...
            par.sigma*par.em*((Temp(i+1,j)+273.15)^4-(par.EL(j).Ta+273.15)^4));             %heat loss to surrounding in the electrolyzer, [watts]
        Qgen(i,j) = par.EL(j).nc*(U(i,j)-par.EL(j).Utn)*I(i,j);                         %heat generated in the electrolyzer, [watts]
        Qlosslye(i,j) = qlye(i,j)*par.Const.CpLye*(Tin(i)-Temp(i+1,j));                 %heat taken out by the lye from the electrolyzer, [watts]
        V_H2(i,j) = nH2elout(i,j)*0.0224136*3600;%hydrogen production rate from individual electrolyzer, [Nm3/h]
        Ps(i,j) = P(i,j)/(1000*V_H2(i,j));%Specific electricity consumption, [kWh/Nm3]
    end
    
    nH2inSto(i) = (nH2in(i)*0.0224136*3600); %net hydrogen flowrate into the storage at all timestamps, [Nm3/h]
    P_net(i)=sum(P(i,:));
    
    if rem(i,100)==0
        disp(mai)
    end
    
    %PI controller in discrete form
    if i >= tstep
        ek(i) = Tin(i)-T_ini(i);
        Ts = 0.1;
        tauC = 2;
        %tau1 = 100;
        Kc = tau1/tauC;
        tauI = min(tauC,4*tau1);
        
        q_cw(i+1) = q_cw(i) + Kc*Ts*ek(i)/tauI + Kc*ek(i);
    end
end

%% Plotting the results
figure()
plot(I)
xlabel('Time, s')
ylabel('Current, A')
legend('El 1','El 2', 'El 3')
grid on

figure()
plot(U)
xlabel('Time, s')
ylabel('Cell voltage, V/cell')
legend('El 1','El 2', 'El 3')
grid on

% figure()
% subplot(1,4,1)
% plot(P(:,1),'b')
% xlabel('Time, s')
% ylabel('P_1, Watts')
% grid on
% subplot(1,4,2)
% plot(P(:,2),'k')
% xlabel('Time, s')
% ylabel('P_2, Watts')
% grid on
% subplot(1,4,3)
% plot(P(:,3),'r')
% xlabel('Time, s')
% ylabel('P_3, Watts')
% grid on
% subplot(1,4,4)
% plot(P_net)
% xlabel('Time, s')
% ylabel('Power, Watts')
% grid on

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
plot(Tin)
subplot(2,1,2)
plot(Tw_out)

figure()
subplot(2,1,1)
plot(PstoH2)
xlabel('Time, s')
ylabel('H_2 Storage pressure, bar')
grid on
subplot(2,1,2)
plot(PcompH2)
xlabel('Time, s')
ylabel('H_2 compressor power, watts')
grid on

figure()
subplot(2,1,1)
plot(Temp(2:end,1))
hold on
plot(Temp(2:end,2))
hold on
plot(Temp(2:end,3))
xlabel('Time, s')
ylabel('T_k, C')
grid on
subplot(2,1,2)
plot(Tout)
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
plot(Tin,'k')
hold on
plot(T_ini,'--')
xlabel('Time, s')
ylabel('T_i_n, C')
ylim([63.5 65.5])
grid on



%% Creating the data file
%save('data_zH2step300_1hr')
%load data_q1step3000_MVQcool_12hr