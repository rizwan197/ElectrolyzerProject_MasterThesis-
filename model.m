function[xDiff, xAlg, input, eqnAlg, eqnDiff] = model(N)
%This function file contains the mathematical model for the system of electrolyzers
%nEl = sequence of the electrolyzer

%Electrolyzer model is a system of following ODE and nonlinear algebraic eqns
%Nonlinear algebraic equation are:
%1N)UI*nc-Power = 0;
%2N)U - (((r1+r2*T)*I)/A) - s*log10(((t1+(t2/T)+(t3/T^2))*I/A)+1) - Urev = 0;
%3N)U*nc - V = 0; U=cell voltage; V=electrolyzer voltage
%4N)Feff - ((.1*I/A)^2)/(f1+((.1*I/A)^2))*f2;
%5N)nH2el - Feff*nc*I/(ze*FC);
%6N)qH2Oloss - nH2*MwtH2O = 0;
%7)nH2_el_net - sum(nH2el) = 0;
%8)nH2_out_net - kvlvH2*VdispH2*sqrt(PstoH2-PoutH2) = 0;
%9)nO2_el_net - nH2elnet/2 = 0;
%10)nO2_out_net - kvlvO2*VdispO2*sqrt(PstoO2-PoutO2) = 0;
%11)T_el_out - ((sum(qlye*T)*CpLye - sum(qloss*T)*Cp + sum(qloss)*(Cp-CpLye)*Tref)/((sum(qlye)-sum(qloss))*Cp) = 0;

%ODE eqautions are:
%1N)dT/dt = (q_lye_k(nEl)*CpLye*(T_El_in-T_k(nEl)) + nc*(u_k(nEl)-cUtn)*i_k(nEl) - ...
%        A_surf*i_k(nEl)* nc/1000*(hc*(T_k(nEl)-par.EL(nEl).Ta) + ...
%        sigma*em*((T_k(nEl)+273.15)^4-(par.EL(nEl).Ta+273.15)^4)))/(CtS*i_k(nEl));
%2)dPstoH2/dt = (TstoH2*Rg/VstoH2)*(nH2-nH2out);
%3)dPstoO2/dt = (TstoO2*Rg/VstoO2)*(nO2-nO2out);
%4)dM_bt/dt = (qlye-qloss) + qH2O - qlye;
%5)dT_El_in/dt = ((netqlye*(T_El_out-T_El_in))/(1000*rho*Vh)) ...
% - (UA*deltaT_LMTD/(1000*rho*CpLye*Vh));
%6)dT_cw_out/dt = ((q_cw*(Tw_in-T_cw_out))/(1000*rho*Vc))...
% + (UA*deltaT_LMTD/(1000*rho*Cp*Vc));

%Inputs for the model are:
%1N)V_el, voltage across each electrolyzer
%2N)qlye_k, lye flowrate through each electrolyzer
%3)q_cw, cooling watwer flowrate 
%4)zH2,outlet valve opening of hydrogen storage tank 
%5)zO2, outlet valve opening of the oxygen storage tank
%6)qH2O, flow rate of water added to buffer tank

%parameters for the model are in parElectrolyzer.m file


%% Load parameters
par = parElectrolyzer(N);

%% Load CasADi
import casadi.*

%% Define symbolic variables
x = SX.sym('x',7*par.N+10);              %symbolic variables for cell voltage, current and electrolyzer voltage (V)
eqnAlg = SX.zeros(6*par.N+5,1);
eqnDiff = SX.zeros(par.N+5,1);

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
T_El_in=x(7*par.N+9);                   %temp of inlet lye stream coming into the electrolyzer
T_cw_out=x(7*par.N+10);


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

inp = SX.sym('inp',2*par.N+4);

U_El_k=[];
q_lye_k=[];
for nEl = 1:par.N
    U_El_k = [U_El_k inp(nEl)];
    q_lye_k = [q_lye_k inp(par.N+nEl)];
end

q_cw = inp(2*par.N+1);
zH2 = inp(2*par.N+2);
zO2 = inp(2*par.N+3);
q_H2O = inp(2*par.N+4);
 

%% Model equations
for nEl = 1:par.N
    
    %relation for Urev with T from LeRoy eqn. 58
    par.EL(nEl).Urev = 1.5184 - 1.5421e-3*(273+T_k(nEl)) + 9.523e-5*(273+T_k(nEl))*log((273+T_k(nEl))) + ...
        9.84e-8*(273+T_k(nEl))^2;
    
    %model equations
    eqnAlg(nEl) = u_k(nEl)*i_k(nEl)*par.EL(nEl).nc - P_k(nEl);                        %power = nc*UI
    eqnAlg(par.N+nEl) = u_k(nEl) - (par.U(nEl).r1+par.U(nEl).r2*T_k(nEl))*i_k(nEl)/par.EL(nEl).A - par.U(nEl).s*log10(((par.U(nEl).t1+par.U(nEl).t2/T_k(nEl)+...
        par.U(nEl).t3/(T_k(nEl)^2))*i_k(nEl)/par.EL(nEl).A)+1) - par.EL(nEl).Urev;  %U-I relationship
    eqnAlg(2*par.N+nEl) = u_k(nEl)*par.EL(nEl).nc - U_El_k(nEl);                    %U(j).nc(j)=V
    eqnAlg(3*par.N+nEl) = Feff_k(nEl) - ((.1*i_k(nEl)/par.EL(nEl).A)^2)/(par.U(nEl).f1+((.1*i_k(nEl)/par.EL(nEl).A)^2))*par.U(nEl).f2;     %faraday efficiency
    eqnAlg(4*par.N+nEl) = nH2_k(nEl) - Feff_k(nEl)*par.EL(nEl).nc*i_k(nEl)/(par.Const.ze*par.Const.FC);%nH2, H2 production rate from individual electrolyzer
    eqnAlg(5*par.N+nEl) = qH2Oloss_k(nEl) - nH2_k(nEl)*par.Const.Mwt;  %flowrate of water lost, [g/s]  
end

sum_H2net = SX.zeros(1,1);
netqlye = SX.zeros(1,1);
netqloss = SX.zeros(1,1);
qlyeT = SX.zeros(1,1);
qlossT = SX.zeros(1,1);

for nEl = 1:par.N
    sum_H2net = sum_H2net + nH2_k(nEl);       %sum of hydrogen from all individual electrolyzers
    netqlye = netqlye + q_lye_k(nEl);             %sum the lye flowing into all individual electrolyzers
    netqloss = netqloss + qH2Oloss_k(nEl);           %sum of water lost from all individual electrolyzers
    qlyeT = qlyeT + q_lye_k(nEl)*T_k(nEl);              %calculate term qlye(k).T(k)
    qlossT = qlossT + qH2Oloss_k(nEl)*T_k(nEl);            %calculate term qH2Oloss(k).T(k)
end


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
eqnDiff(par.N+3) = (netqlye-netqloss) + q_H2O - netqlye;                                     %differential eqn for mass in the buffer tank, [grams i.e. pho*V]

%dynamic thermal balance for the heat exchanger, T_in and Tw_out are differential variables

deltaT1 = T_El_in - par.Tw_in;%difference in temperatures of hot and cold streams at inlet
deltaT2 = T_El_out - T_cw_out;%difference in temperatures of hot and cold streams at outlet
deltaT_LMTD = (deltaT1-deltaT2)/log(deltaT1/deltaT2);

eqnDiff(par.N+4) = ((netqlye*(T_El_out-T_El_in))/(1000*par.Const.rho*par.Const.Vh)) - (par.Hex.UA*deltaT_LMTD/(1000*par.Const.rho*par.Const.CpLye*par.Const.Vh));%differential eqn for the hot stream exit temp from heat exchanger
eqnDiff(par.N+5) = ((q_cw*(par.Tw_in-T_cw_out))/(1000*par.Const.rho*par.Const.Vc)) + (par.Hex.UA*deltaT_LMTD/(1000*par.Const.rho*par.Const.Cp*par.Const.Vc));%differential eqn for the cold stream exit temperature from heat exchanger

xDiff = x(6*par.N+6:end,1);
xAlg = x(1:6*par.N+5,1);
input = inp;
end
