function[xDiff, xAlg, input, eqnAlg, eqnDiff, F] = modelPI(N,u0,Kc,tauI,PC,T2C,T3C)
%This function file contains the mathematical model for the system of
%electrolyzers + PI controller
%This plant model includes extra states for the integrated error of the CVs 
%nEl = sequence of the electrolyzer
%Electrolyzer model is a system of following ODE and nonlinear algebraic eqns
%Nonlinear algebraic equation are:
%1N)UI*nc-Power = 0;
%2N)U - (((r1+r2*T)*I)/A) - s*log10(((t1+(t2/T)+(t3/T^2))*I/A)+1) - Urev = 0;
%3N)U*nc - V = 0; U=cell voltage; V=electrolyzer voltage
%4N)Feff - ((.1*I/A)^2)/(f1+((.1*I/A)^2))*f2;
%5N)nH2el - Feff*nc*I/(ze*FC);
%6N)qH2Oloss - nH2*MwtH2O = 0;
%6N+1)nH2_el_net - sum(nH2el) = 0;
%6N+2)nH2_out_net - kvlvH2*VdispH2*sqrt(PstoH2-PoutH2) = 0;
%6N+3)nO2_el_net - nH2elnet/2 = 0;
%6N+4)nO2_out_net - kvlvO2*VdispO2*sqrt(PstoO2-PoutO2) = 0;
%6N+5)T_el_out - ((sum(qlye*T)*CpLye - sum(qloss*T)*Cp + sum(qloss)*(Cp-CpLye)*Tref)/((sum(qlye)-sum(qloss))*Cp) = 0;
%ODE eqautions are:
%1N)dT/dt = (q_lye_k(nEl)*CpLye*(T_El_in-T_k(nEl)) + nc*(u_k(nEl)-cUtn)*i_k(nEl) - ...
%        A_surf*(hc*(T_k(nEl)-par.EL(nEl).Ta) + ...
%        sigma*em*((T_k(nEl)+273.15)^4-(par.EL(nEl).Ta+273.15)^4)))/(CtS*Pnom);
%N+1)dPstoH2/dt = (TstoH2*Rg/VstoH2)*(nH2-nH2out);
%N+2)dPstoO2/dt = (TstoO2*Rg/VstoO2)*(nO2-nO2out);
%N+3)dM_bt/dt = (qlye-qloss) + qH2O - qlye;
%N+4)dT_bt_out/dt = (netqout*par.Const.CpLye*(T_El_out-T_bt_out) + ...
%     q_H2O*(par.Const.Cp*T_bt_out - par.Const.CpLye*T_bt_out) ...
%     - (netqout*par.Const.CpLye+q_H2O*par.Const.Cp-netqlye*par.Const.CpLye)*par.Const.Tref)...
%     /(Mass_Bt*par.Const.CpLye);
%N+5)dT_El_in/dt = ((netqlye*(T_bt_out-T_El_in))/(1000*rho*Vh)) ...
%                - (UA*deltaT_LMTD/(1000*rho*CpLye*Vh));
%N+6)dT_cw_out/dt = ((q_cw*(Tw_in-T_cw_out))/(1000*rho*Vc))...
%                 + (UA*deltaT_LMTD/(1000*rho*Cp*Vc));
%Inputs for the model are:
%1N)V_el, voltage across each electrolyzer
%2N)qlye_k, lye flowrate through each electrolyzer
%3)q_cw, cooling water flowrate 
%4)zH2,outlet valve opening of hydrogen storage tank 
%5)zO2, outlet valve opening of the oxygen storage tank
%6)qH2O, flow rate of water added to buffer tank
%parameters for the model are in parElectrolyzer.m file


%% Load parameters
par = parElectrolyzer(N);

%% Load CasADi
import casadi.*

%% Define symbolic variables
x = SX.sym('x',7*par.N+16);              %symbolic variables for cell voltage, current and electrolyzer voltage (V)
eqnAlg = SX.zeros(6*par.N+5,1);
eqnDiff = SX.zeros(par.N+11,1);

%variables for algebriac eqns.
u_k=[];
i_k=[];
P_k=[];
Feff_k=[];
nH2_k=[];
qH2Oloss_k=[];
nH2El_net=x(6*par.N+1);
nH2out_net=x(6*par.N+2);
nO2El_net=x(6*par.N+3);
nO2out_net=x(6*par.N+4);
T_El_out=x(6*par.N+5);                  %temp after mixing of the exiting liquid streams from all electrolyzers

%variables for differential equations 
T_k=[];
Psto_H2=x(7*par.N+6);
Psto_O2=x(7*par.N+7);
Mass_Bt=x(7*par.N+8);
T_bt_out=x(7*par.N+9);
T_El_in=x(7*par.N+10);                   %temp of inlet lye stream coming into the electrolyzer
T_cw_out=x(7*par.N+11);

for nEl = 1:par.N
    u_k = [u_k x(nEl)];                         %cell voltage of the electrolyzer
    i_k = [i_k x(par.N+nEl)];                   %current in the electrolyzer
    P_k = [P_k x(2*par.N+nEl)];                 %power of the individual electrolyzer
    Feff_k = [Feff_k x(3*par.N+nEl)];           %faraday efficiency of each electrolyzer
    nH2_k = [nH2_k x(4*par.N+nEl)];             %hydrogen produced form each individual electrolyzer
    qH2Oloss_k = [qH2Oloss_k x(5*par.N+nEl)];   %water loss during electrolysis in kth electrolyzer
    T_k = [T_k x(6*par.N+5+nEl)];               %temperature of the individual electrolyzer
end 

%% Define inputs for the simulation (MVs for the dynamic simulation)

inp = SX.sym('inp',par.N+5);
q_lye_k=[];

Pnetset = inp(1);
q_lye_k(1) = inp(2);
T2set = inp(3);
T3set = inp(4);

q_cw = inp(par.N+2);
pstoH2set = inp(par.N+3);
pstoO2set = inp(par.N+4);
Mass_Btset = inp(par.N+5);

%% PI controller for states with zero steady state effect

%PI controller for pressure in the hydrogen storage tank
eint_pstoH2 = x(7*par.N+12);            %eint in pstoH2
e_pstoH2 = Psto_H2 - pstoH2set;         %error in pstoH2
Kc_pstoH2PI = Kc(1);                    %controller gain
taui_pstoH2PI = tauI(1);                %integral time constant
zH2ini = u0(1);                         %initial value of MV i.e. zH2
zH2 = min(1,max(0,PIcontroller(zH2ini,Kc_pstoH2PI,taui_pstoH2PI,e_pstoH2,eint_pstoH2)));%0<=zO2<=1

%PI controller for pressure in the oxygen storage tank
eint_pstoO2 = x(7*par.N+13);            %eint in mass of the buffer tank
e_pstoO2 = Psto_O2 - pstoO2set;         %error
Kc_pstoO2PI = Kc(2);                    %controller gain
taui_pstoO2PI = tauI(2);                %integral time constant
zO2ini = u0(2);                         %initial value of MV i.e. zO2
zO2 = min(1,max(0,PIcontroller(zO2ini,Kc_pstoO2PI,taui_pstoO2PI,e_pstoO2,eint_pstoO2)));%0<=zO2<=1

%PI controller for mass in the buffer tank
eint_Mbt = x(7*par.N+14);               %eint in mass of the buffer tank
eMass_Bt = Mass_Bt - Mass_Btset;        %error
Kc_MassBtPI = Kc(3);                    %controller gain
taui_MassBtPI = tauI(3);                %integral time constant
q_H2Oini = u0(3);                       %initial value of MV i.e. qH2O
q_H2O = max(0,PIcontroller(q_H2Oini,Kc_MassBtPI,taui_MassBtPI,eMass_Bt,eint_Mbt));%q_H2O>=0

%% Controller for rest of the states

%proportional controller for Vel-Pnet
Pcons = SX.zeros(1,1);
for nEl = 1:par.N
    Pcons = Pcons+P_k(nEl);
end
PC.e = Pcons-Pnetset;
U_El_k =  PC.u0 - (PC.Kc)*PC.e;

%PI controller for the temperature of the electrolyzer 2
eint_T2 = x(7*par.N+15);     %eint T2
eT2 = T_k(2) - T2set;        %error
q_lye_k(2) = min(10,max(0.5,PIcontroller(T2C.u0,T2C.Kc,T2C.tauI,eT2,eint_T2)));%0.5<=q_lye_k<=10;

%PI controller for the temperature of the electrolyzer 3
eint_T3 = x(7*par.N+16);     %eint T3
eT3 = T_k(3) - T3set;        %error
q_lye_k(3) = min(10,max(0.5,PIcontroller(T3C.u0,T3C.Kc,T3C.tauI,eT3,eint_T3)));%0.5<=q_lye_k<=10;

%% Model equations
for nEl = 1:par.N
    
    %relation for Urev with T from LeRoy eqn. 58
    par.EL(nEl).Urev = 1.5184 - 1.5421e-3*(273+T_k(nEl)) + 9.523e-5*...
        (273+T_k(nEl))*log((273+T_k(nEl))) + 9.84e-8*(273+T_k(nEl))^2;
    
    %model equations
    eqnAlg(nEl) = u_k(nEl)*i_k(nEl)*par.EL(nEl).nc - P_k(nEl);                      %power = nc*UI
    eqnAlg(par.N+nEl) = u_k(nEl) - (par.U(nEl).r1+par.U(nEl).r2*...
        T_k(nEl))*i_k(nEl)/par.EL(nEl).A - par.U(nEl).s*...
        log10(((par.U(nEl).t1+par.U(nEl).t2/T_k(nEl)+...
        par.U(nEl).t3/(T_k(nEl)^2))*i_k(nEl)/par.EL(nEl).A)+1)-...
        par.EL(nEl).Urev;                                                           %U-I relationship
    eqnAlg(2*par.N+nEl) = u_k(nEl)*par.EL(nEl).nc - U_El_k;                    %U(j).nc(j)=V
    eqnAlg(3*par.N+nEl) = Feff_k(nEl) - ((.1*i_k(nEl)/par.EL(nEl).A)^2)/...
        (par.U(nEl).f1+((.1*i_k(nEl)/par.EL(nEl).A)^2))*par.U(nEl).f2;              %faraday efficiency
    eqnAlg(4*par.N+nEl) = nH2_k(nEl) - Feff_k(nEl)*par.EL(nEl).nc*...
        i_k(nEl)/(par.Const.ze*par.Const.FC);                                       %nH2, H2 production rate from individual electrolyzer
    eqnAlg(5*par.N+nEl) = qH2Oloss_k(nEl) - nH2_k(nEl)*par.Const.Mwt;               %flowrate of water lost, [g/s]  
end

sum_H2net = SX.zeros(1,1);
netqlye = SX.zeros(1,1);
netqloss = SX.zeros(1,1);
qlyeT = SX.zeros(1,1);
qlossT = SX.zeros(1,1);

for nEl = 1:par.N
    sum_H2net = sum_H2net + nH2_k(nEl);         %sum of hydrogen from all individual electrolyzers
    netqlye = netqlye + q_lye_k(nEl);           %sum the lye flowing into all individual electrolyzers
    netqloss = netqloss + qH2Oloss_k(nEl);      %sum of water lost from all individual electrolyzers
    qlyeT = qlyeT + q_lye_k(nEl)*T_k(nEl);      %calculate term qlye(k).T(k)
    qlossT = qlossT + qH2Oloss_k(nEl)*T_k(nEl); %calculate term qH2Oloss(k).T(k)
end
netqout = netqlye-netqloss;

eqnAlg(6*par.N+1) = nH2El_net - sum_H2net;                                              %algebraic eqn for net hydrogen flowrate from all the electrolyzers
eqnAlg(6*par.N+2) = nH2out_net - par.kvalveH2*zH2*sqrt(Psto_H2-par.Storage.PoutH2);     %algebraic eqn for net hydrogen flowrate from the storage tank
eqnAlg(6*par.N+3) = nO2El_net - nH2El_net/2;                                            %algebraic eqn for net oxygen flowrate from all the electrolyzers
eqnAlg(6*par.N+4) = nO2out_net - par.kvalveO2*zO2*sqrt(Psto_O2-par.Storage.PoutO2);     %algebraic eqn for net oxygen flowrate from the storage tank
eqnAlg(6*par.N+5) = T_El_out - (((qlyeT*par.Const.CpLye) - (qlossT*par.Const.Cp) + netqloss*(par.Const.Cp-par.Const.CpLye)*par.Const.Tref)/...
    ((netqlye-netqloss)*par.Const.CpLye));                                          %calculation of Tout

for nEl = 1:par.N
    eqnDiff(nEl) = (q_lye_k(nEl)*par.Const.CpLye*(T_El_in-T_k(nEl)) + par.EL(nEl).nc*(u_k(nEl)-par.EL(nEl).Utn)*i_k(nEl) - ...
        par.TherMo(nEl).A_El*(par.TherMo(nEl).hc*(T_k(nEl)-par.EL(nEl).Ta) + ...
        par.sigma*par.em*((T_k(nEl)+273.15)^4-(par.EL(nEl).Ta+273.15)^4)))/(par.TherMo(nEl).Ct*1000);%differential eqn for the electrolyzer temperature
end

eqnDiff(par.N+1) = (par.Storage.TstoH2*par.Storage.Rg/par.Storage.VstoH2)*(nH2El_net-nH2out_net);  %differential eqn for hydrogen storage pressure
eqnDiff(par.N+2) = (par.Storage.TstoO2*par.Storage.Rg/par.Storage.VstoO2)*(nO2El_net-nO2out_net);  %differential eqn for oxygen storage pressure
eqnDiff(par.N+3) = netqout + q_H2O - netqlye;    %differential eqn for mass of liquid in the buffer tank, [grams i.e. pho*V]
%the level remains same at steady state but starts to change with the change
%in net hydrogen production, since the netqloss changes whereas q_H2O is a
%MV (constant parameter for intergration over time for dynamic states).
eqnDiff(par.N+4) = (netqout*par.Const.CpLye*(T_El_out-T_bt_out) + ...
    q_H2O*(par.Const.Cp*T_bt_out - par.Const.CpLye*T_bt_out) ...
    - (netqout*par.Const.CpLye+q_H2O*par.Const.Cp-netqlye*par.Const.CpLye)*par.Const.Tref)/(Mass_Bt*par.Const.CpLye);
%assuming T_H2O = T_bt_out

%dynamic thermal balance for the heat exchanger, T_in and Tw_out are differential variables
deltaT1 = T_El_in - par.Tw_in;%difference in temperatures of hot and cold streams at inlet
deltaT2 = T_bt_out - T_cw_out;%difference in temperatures of hot and cold streams at outlet
deltaT_LMTD = (deltaT1-deltaT2)/log(deltaT1/deltaT2);

eqnDiff(par.N+5) = ((netqlye*(T_bt_out-T_El_in))/(1000*par.Const.rhoLye*par.Const.Vh)) - ...
    (par.Hex.UA*deltaT_LMTD/(1000*par.Const.rhoLye*par.Const.CpLye*par.Const.Vh));%differential eqn for the hot stream exit temp from heat exchanger
eqnDiff(par.N+6) = ((q_cw*(par.Tw_in-T_cw_out))/(1000*par.Const.rho*par.Const.Vc)) + ...
    (par.Hex.UA*deltaT_LMTD/(1000*par.Const.rho*par.Const.Cp*par.Const.Vc));%differential eqn for the cold stream exit temperature from heat exchanger
eqnDiff(par.N+7) = e_pstoH2; %differential equation d(eint_pstoH2)/dt = e_pstoH2
eqnDiff(par.N+8) = e_pstoO2; %differential equation d(eint_pstoO2)/dt = e_pstoO2
eqnDiff(par.N+9) = eMass_Bt; %differential equation d(eintMass_Bt)/dt = eMass_Bt
eqnDiff(par.N+10) = eint_T2;
eqnDiff(par.N+11) = eint_T3;

xDiff = x(6*par.N+6:end,1);
xAlg = x(1:6*par.N+5,1);
input = inp;

dae = struct('x',xDiff,'z',xAlg,'p',input,'ode',eqnDiff,'alg',eqnAlg);
F = integrator('F', 'idas', dae);

end