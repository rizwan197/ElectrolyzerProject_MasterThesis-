function[xDiff, xAlg, input, eqnAlg, eqnDiff, F] = model(N)
%This function file contains the mathematical model for the system of electrolyzers
%nEl = sequence of the electrolyzer

%Electrolyzer model is a system of following ODE and nonlinear algebraic eqns
%Nonlinear algebraic equation are:
%1N)UI*nc-Power = 0;
%2N)U - (((r1+r2*T)*I)/A) - s*log10(((t1+(t2/T)+(t3/T^2))*I/A)+1) - Urev = 0;
%3N)U*nc - V = 0; U=cell voltage; V=electrolyzer voltage
%4N)Feff - ((.1*I/A)^2)/(f1+((.1*I/A)^2))*f2;
%5N)nH2el - Feff*nc*I/(ze*FC);

%ODE eqautions are:
%1N)dT/dt = (q_lye_k(nEl)*CpLye*(T_El_in-T_k(nEl)) + nc*(u_k(nEl)-cUtn)*i_k(nEl) - ...
%        A_surf*(hc*(T_k(nEl)-par.EL(nEl).Ta) + ...
%        sigma*em*((T_k(nEl)+273.15)^4-(par.EL(nEl).Ta+273.15)^4)))/(CtS*Pnom);

%Inputs for the model are:
%1N)qlye_k, lye flowrate through each electrolyzer

%parameters for the model are in parElectrolyzer.m file


%% Load parameters
par = parElectrolyzer(N);

%% Load CasADi
import casadi.*

%% Define symbolic variables
x = SX.sym('x',6*par.N);              %symbolic variables for cell voltage, current and electrolyzer voltage (V)
eqnAlg = SX.zeros(5*par.N,1);
eqnDiff = SX.zeros(par.N,1);

%variables for algebriac eqns.
u_k=[];
i_k=[];
P_k=[];
Feff_k=[];
nH2_k=[];

%variables for differential equations 
T_k=[];
T_El_in = 65;

for nEl = 1:par.N
    %algerbriac variables
    u_k = [u_k x(nEl)];                         %cell voltage of the electrolyzer
    i_k = [i_k x(par.N+nEl)];                   %current in the electrolyzer
    P_k = [P_k x(2*par.N+nEl)];                 %power of the individual electrolyzer
    Feff_k = [Feff_k x(3*par.N+nEl)];           %faraday efficiency of each electrolyzer
    nH2_k = [nH2_k x(4*par.N+nEl)];             %hydrogen produced form each individual electrolyzer
    %differential variables
    T_k = [T_k x(5*par.N+nEl)];               %temperature of the individual electrolyzer
end 

%% Define inputs for the simulation (MVs for the dynamic simulation)

inp = SX.sym('inp',par.N);

q_lye_k=[];
for nEl = 1:par.N
    q_lye_k = [q_lye_k inp(nEl)];
end
 

%% Model equations
for nEl = 1:par.N
    
    %relation for Urev with T from LeRoy eqn. 58
    par.EL(nEl).Urev = 1.5184 - 1.5421e-3*(273+T_k(nEl)) + 9.523e-5*(273+T_k(nEl))*log((273+T_k(nEl))) + ...
        9.84e-8*(273+T_k(nEl))^2;
    
    %model equations
    eqnAlg(nEl) = u_k(nEl)*i_k(nEl)*par.EL(nEl).nc - P_k(nEl);                        %power = nc*UI
    eqnAlg(par.N+nEl) = u_k(nEl) - (par.U(nEl).r1+par.U(nEl).r2*T_k(nEl))*i_k(nEl)/par.EL(nEl).A - par.U(nEl).s*log10(((par.U(nEl).t1+par.U(nEl).t2/T_k(nEl)+...
        par.U(nEl).t3/(T_k(nEl)^2))*i_k(nEl)/par.EL(nEl).A)+1) - par.EL(nEl).Urev;  %U-I relationship
    
    eqnAlg(2*par.N+nEl) = Feff_k(nEl) - ((.1*i_k(nEl)/par.EL(nEl).A)^2)/(par.U(nEl).f1+((.1*i_k(nEl)/par.EL(nEl).A)^2))*par.U(nEl).f2;     %faraday efficiency
    eqnAlg(3*par.N+nEl) = nH2_k(nEl) - Feff_k(nEl)*par.EL(nEl).nc*i_k(nEl)/(par.Const.ze*par.Const.FC);%nH2, H2 production rate from individual electrolyzer
    
end


for nEl = 1:par.N
    eqnDiff(nEl) = (q_lye_k(nEl)*par.Const.CpLye*(T_El_in-T_k(nEl)) + par.EL(nEl).nc*(u_k(nEl)-par.EL(nEl).Utn)*i_k(nEl) - ...
        par.TherMo(nEl).A_El*(par.TherMo(nEl).hc*(T_k(nEl)-par.EL(nEl).Ta) + ...
        par.sigma*par.em*((T_k(nEl)+273.15)^4-(par.EL(nEl).Ta+273.15)^4)))/(par.TherMo(nEl).Ct*1000);%differential eqn for the electrolyzer temperature
end


xDiff = x(5*par.N+1:end,1);
xAlg = x(1:5*par.N,1);
input = inp;

dae = struct('x',xDiff,'z',xAlg,'p',input,'ode',eqnDiff,'alg',eqnAlg);
F = integrator('F', 'idas', dae);

end