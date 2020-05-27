function[xDiff, xAlg, input, eqnAlg, eqnDiff, F] = modelEdit(N)
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
%7N)T_El_out_k(nEl) - ((q_lye_k(nEl)*T_k(nEl)*par.Const.CpLye - qH2Oloss_k(nEl)*T_k(nEl)*par.Const.Cp + qH2Oloss_k(nEl)*(par.Const.CpLye-par.Const.Cp)*...
       % par.Const.Tref)/(q_lye_k(nEl)*par.Const.CpLye - qH2Oloss_k(nEl)*par.Const.CpLye)) = 0 ;

%ODE eqautions are:
%1N)dT/dt = (q_lye_k(nEl)*CpLye*(T_El_in-T_k(nEl)) + nc*(u_k(nEl)-cUtn)*i_k(nEl) - ...
%        A_surf*(hc*(T_k(nEl)-par.EL(nEl).Ta) + ...
%        sigma*em*((T_k(nEl)+273.15)^4-(par.EL(nEl).Ta+273.15)^4)))/(CtS*Pnom);
%2N)dM_bt/dt = (qlye-qloss) + qH2O - qlye;
%3N)dT_bt_out/dt = (netqout*par.Const.CpLye*(T_El_out-T_bt_out) + ...
%     q_H2O*(par.Const.Cp*T_bt_out - par.Const.CpLye*T_bt_out) ...
%     - (netqout*par.Const.CpLye+q_H2O*par.Const.Cp-netqlye*par.Const.CpLye)*par.Const.Tref)...
%     /(Mass_Bt*par.Const.CpLye);
%4N)dT_El_in/dt = ((netqlye*(T_bt_out-T_El_in))/(1000*rho*Vh)) ...
%                - (UA*deltaT_LMTD/(1000*rho*CpLye*Vh));
%5N)dT_cw_out/dt = ((q_cw*(Tw_in-T_cw_out))/(1000*rho*Vc))...
%                 + (UA*deltaT_LMTD/(1000*rho*Cp*Vc));

%Inputs for the model are:
%1N)V_el, voltage across each electrolyzer
%2N)qlye_k, lye flowrate through each electrolyzer
%3N)q_cw, cooling water flowrate
%4N)qH2O, flow rate of water added to buffer tank

%parameters for the model are in parElectrolyzer.m file


%% Load parameters
par = parElectrolyzer(N);

%% Load CasADi
import casadi.*

%% Define symbolic variables
x = SX.sym('x',12*par.N);              %symbolic variables for cell voltage, current and electrolyzer voltage (V)
eqnAlg = SX.zeros(9*par.N,1);
eqnDiff = SX.zeros(3*par.N,1);

%variables for algebriac eqns.
u_k=[];
i_k=[];
P_k=[];
Feff_k=[];
nH2_k=[];
qH2Oloss_k=[];
T_bt_in_k=[];
T_El_in_k=[];                   %temp of inlet lye stream coming into the electrolyzer
T_cw_out_k=[];

%variables for differential equations
T_k=[];
Mass_Bt_k=[];
T_bt_out_k=[];


for nEl = 1:par.N
    %algebriac variables
    u_k = [u_k x(nEl)];                         %cell voltage of the electrolyzer
    i_k = [i_k x(par.N+nEl)];                   %current in the electrolyzer
    P_k = [P_k x(2*par.N+nEl)];                 %power of the individual electrolyzer
    Feff_k = [Feff_k x(3*par.N+nEl)];           %faraday efficiency of each electrolyzer
    nH2_k = [nH2_k x(4*par.N+nEl)];             %hydrogen produced form each individual electrolyzer
    qH2Oloss_k = [qH2Oloss_k x(5*par.N+nEl)];   %water loss during electrolysis in kth electrolyzer
    T_bt_in_k=[T_bt_in_k x(6*par.N+nEl)];       %temperature of lye stream entering into the buffer tank
    T_El_in_k=[T_El_in_k x(7*par.N+nEl)];       %temp of inlet lye stream coming into the electrolyzer
    T_cw_out_k=[T_cw_out_k x(8*par.N+nEl)];
    
    %differential variables
    T_k = [T_k x(9*par.N+nEl)];                 %temperature of the individual electrolyzer
    Mass_Bt_k=[Mass_Bt_k x(10*par.N+nEl)];
    T_bt_out_k=[T_bt_out_k x(11*par.N+nEl)];    %temperature of lye stream leaving the buffer tank
    
end

%% Define inputs for the simulation (MVs for the dynamic simulation)

inp = SX.sym('inp',4*par.N);

U_El_k=[];
q_lye_k=[];
Q_cool_k=[];%[J/s]
q_H2O_k=[];

for nEl = 1:par.N
    U_El_k = [U_El_k inp(nEl)];
    q_lye_k = [q_lye_k inp(par.N+nEl)];
    Q_cool_k = [Q_cool_k inp(2*par.N+nEl)];
    q_H2O_k = [q_H2O_k inp(3*par.N+nEl)];
end

%% Model equations
for nEl = 1:par.N
    
    %relation for Urev with T from LeRoy eqn. 58
    par.EL(nEl).Urev = 1.5184 - 1.5421e-3*(273+T_k(nEl)) + 9.523e-5*(273+T_k(nEl))*log((273+T_k(nEl))) + ...
        9.84e-8*(273+T_k(nEl))^2;
    
    %model equations
    eqnAlg(nEl) = u_k(nEl)*i_k(nEl)*par.EL(nEl).nc - P_k(nEl);                          %power = nc*UI
    eqnAlg(par.N+nEl) = u_k(nEl) - (par.U(nEl).r1+par.U(nEl).r2*T_k(nEl))*i_k(nEl)/par.EL(nEl).A - par.U(nEl).s*log10(((par.U(nEl).t1+par.U(nEl).t2/T_k(nEl)+...
        par.U(nEl).t3/(T_k(nEl)^2))*i_k(nEl)/par.EL(nEl).A)+1) - par.EL(nEl).Urev;      %U-I relationship
    eqnAlg(2*par.N+nEl) = u_k(nEl)*par.EL(nEl).nc - U_El_k(nEl);                        %U(j).nc(j)=V
    eqnAlg(3*par.N+nEl) = Feff_k(nEl) - ((.1*i_k(nEl)/par.EL(nEl).A)^2)/(par.U(nEl).f1+((.1*i_k(nEl)/par.EL(nEl).A)^2))*par.U(nEl).f2;     %faraday efficiency
    eqnAlg(4*par.N+nEl) = nH2_k(nEl) - Feff_k(nEl)*par.EL(nEl).nc*i_k(nEl)/(par.Const.ze*par.Const.FC);%nH2, H2 production rate from individual electrolyzer
    eqnAlg(5*par.N+nEl) = qH2Oloss_k(nEl) - nH2_k(nEl)*par.Const.Mwt;                   %flowrate of water lost, [g/s]
    eqnAlg(6*par.N+nEl) = T_bt_in_k(nEl) - ((q_lye_k(nEl)*T_k(nEl)*par.Const.CpLye - qH2Oloss_k(nEl)*T_k(nEl)*par.Const.Cp + qH2Oloss_k(nEl)*(par.Const.Cp-par.Const.CpLye)*...
        par.Const.Tref)/(q_lye_k(nEl)*par.Const.CpLye - qH2Oloss_k(nEl)*par.Const.CpLye));                  %calculation of T_bt_in_k
    eqnAlg(7*par.N+nEl) = Q_cool_k(nEl) - (q_lye_k(nEl)*par.Const.CpLye*(T_bt_out_k(nEl)-T_El_in_k(nEl)));  %calculation of T_El_in
    eqnAlg(8*par.N+nEl) = Q_cool_k(nEl) + ((8e4/N)*par.Const.Cp*(par.Tw_in_k-T_cw_out_k(nEl)));         %calculation of T_cw_out
end

qout_k = q_lye_k - qH2Oloss_k;

for nEl = 1:par.N
    eqnDiff(nEl) = (q_lye_k(nEl)*par.Const.CpLye*(T_El_in_k(nEl)-T_k(nEl)) + par.EL(nEl).nc*(u_k(nEl)-par.EL(nEl).Utn)*i_k(nEl) - ...
        par.TherMo(nEl).A_El*(par.TherMo(nEl).hc*(T_k(nEl)-par.EL(nEl).Ta) + ...
        par.sigma*par.em*((T_k(nEl)+273.15)^4-(par.EL(nEl).Ta+273.15)^4)))/(par.TherMo(nEl).Ct*1000);%differential eqn for the electrolyzer temperature
    
    eqnDiff(par.N+nEl) = qout_k(nEl) + q_H2O_k(nEl) - q_lye_k(nEl);    %differential eqn for mass of liquid in the buffer tank, [grams i.e. pho*V]
    %the level remains same at steady state but starts to change with the change
    %in net hydrogen production, since the netqloss changes whereas q_H2O_k is a
    %MV (constant parameter for intergration over time for dynamic states).
    
    eqnDiff(2*par.N+nEl) = (qout_k(nEl)*par.Const.CpLye*(T_bt_in_k(nEl)-T_bt_out_k(nEl)) + ...
        q_H2O_k(nEl)*(par.Const.Cp*T_bt_out_k(nEl) - par.Const.CpLye*T_bt_out_k(nEl)) ...
        - (qout_k(nEl)*par.Const.CpLye+q_H2O_k(nEl)*par.Const.Cp-q_lye_k(nEl)*par.Const.CpLye)*par.Const.Tref)/(Mass_Bt_k(nEl)*par.Const.CpLye);
    %assuming T_H2O = T_bt_out
end

% %dynamic thermal balance for the heat exchanger, T_in and Tw_out are differential variables
% deltaT1 = T_El_in_k - par.Tw_in_k*ones(1,par.N);%difference in temperatures of hot and cold streams at inlet
% deltaT2 = T_bt_out_k - T_cw_out_k;%difference in temperatures of hot and cold streams at outlet
% deltaT_LMTD = (deltaT1-deltaT2)./log(deltaT1./deltaT2);
% 
% for nEl=1:par.N  
%     eqnDiff(3*par.N+nEl) = ((q_lye_k(nEl)*(T_bt_out_k(nEl)-T_El_in_k(nEl)))/(1000*par.Const.rhoLye*par.Const.Vhk)) - ...
%         (par.Hex.UAk*deltaT_LMTD(nEl)/(1000*par.Const.rhoLye*par.Const.CpLye*par.Const.Vhk));%differential eqn for the hot stream exit temp from heat exchanger
%     
%     eqnDiff(4*par.N+nEl) = ((q_cw_k(nEl)*(par.Tw_in_k-T_cw_out_k(nEl)))/(1000*par.Const.rho*par.Const.Vck)) + ...
%         (par.Hex.UAk*deltaT_LMTD(nEl)/(1000*par.Const.rho*par.Const.Cp*par.Const.Vck));%differential eqn for the cold stream exit temperature from heat exchanger
% end

xAlg = x(1:9*par.N,1);
xDiff = x(9*par.N+1:end,1);
input = inp;

dae = struct('x',xDiff,'z',xAlg,'p',input,'ode',eqnDiff,'alg',eqnAlg);
F = integrator('F', 'idas', dae);

end