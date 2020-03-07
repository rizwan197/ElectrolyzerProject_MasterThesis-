function [x0,z0] = elsstr(param,par)

%This function file back calculates (numerical solution) the minimum lye flow rate required for
%the operation of electrolyzer plant
%Put a hold at steady state solution in the main file and run this file
%Here we will assume temperature rises in each electrolyzer to 80C at steady state
%Power, qlye and V (electrolyzer voltage) are provided as inputs to solve for steady state

%nEl = sequence of the electrolyzer
%Here we are solving following eqns for each electrolyzer:
%1N)UI*nc-Power = 0;
%2N)U - (((r1+r2*T)*I)/A) - s*log10(((t1+(t2/T)+(t3/T^2))*I/A)+1) - Urev = 0;
%3N)U*nc - V = 0; U=cell voltage; V=electrolyzer voltage
%4N)qlye*CpLye*(T_in-T) + nc*(U-Utn)*I - (1/Rt)*(T-Ta);

%% Load CasADi
%addpath('/Users/mdrizwan/Documents/MATLAB/casadi-osx-matlabR2015a-v3.4.5')
import casadi.*

%% Define symbolic variables
x = SX.sym('x',4*par.N);            %symbolic variables for cell voltage, current and electrolyzer voltage (V)
eqn = SX.zeros(4*par.N,1);
u=[];
i=[];
qlye_k=[];
Ps=[];


for nEl=1:par.N
    u=[u x(nEl)];                   %cell voltage of the electrolyzer
    i=[i x(par.N+nEl)];             %current in the electrolyzer
    qlye_k=[qlye_k x(2*par.N+nEl)];       %temperature of the individual electrolyzer
    Ps=[Ps x(3*par.N+nEl)];         %power of the individual electrolyzer
end

%% Initial conditions for steady state calculations


T_in = param.T0;                                        %steady state temperature of lye into electrolyzer

u0=zeros(1,par.N);
i0=zeros(1,par.N);
Ps0=zeros(1,par.N);
q_lye0=zeros(1,par.N);

T_k0=zeros(1,par.N);

for nEl = 1:par.N
    u0(nEl) = param.u0*21/par.EL(nEl).nc;                       %initial guess for cell voltage
    Ps0(nEl) = param.Ps0;
    q_lye0(nEl) = param.q_lye(nEl);
    i0(nEl) = Ps0(nEl)/(par.EL(nEl).nc*u0(nEl));                %initial guess for current
    T_k0(nEl) = 80;

end

X0 = [u0 i0 q_lye0 Ps0];

%% Solving steady state problem
for nEl = 1:par.N
    eqn(nEl) = u(nEl)*i(nEl)*par.EL(nEl).nc-Ps(nEl);                        %power = nc*UI
    eqn(par.N+nEl) = u(nEl) - (par.U(nEl).r1+par.U(nEl).r2*T_k0(nEl))*i(nEl)/par.EL(nEl).A - par.U(nEl).s*log10(((par.U(nEl).t1+par.U(nEl).t2/T_k0(nEl)+...
        par.U(nEl).t3/(T_k0(nEl)^2))*i(nEl)/par.EL(nEl).A)+1) - par.EL(nEl).Urev;                                              %U-I relationship
    eqn(2*par.N+nEl) = qlye_k(nEl)*par.Const.CpLye*(T_in-T_k0(nEl)) + par.EL(nEl).nc*(u(nEl)-par.EL(nEl).Utn)*i(nEl) - ...
        (1/par.TherMo(nEl).Rt)*(T_k0(nEl)-par.EL(nEl).Ta);                         %temperature of each electrolyzer
end

for nEl = 1:(par.N-1)
    eqn(3*par.N+nEl) = u(nEl)*par.EL(nEl).nc - u(nEl+1)*par.EL(nEl+1).nc;                                %U(j).nc(j)= V
end

sumPs = SX.zeros(1,1);
for nEl = 1:par.N
   sumPs =  sumPs + Ps(nEl);
end

eqn(4*par.N) = (sumPs - sum(Ps0));                                    % V = voltage across last electrolyzer


g = Function('g',{x},{eqn});
G = rootfinder('G','newton',g);
res = full(G(X0));

U=[];
I=[];
Pk=[];
Qlyek=[];
for nEl=1:par.N
    U = [U res(nEl)];
    I = [I res(par.N+nEl)];
    Qlyek = [Qlyek res(2*par.N+nEl)];
    Pk = [Pk res(3*par.N+nEl)];
end
U
I
Qlyek
Pk


%% Residuals check

 for nEl=1:par.N
 eq1(nEl) = I(nEl)*U(nEl)*par.EL(nEl).nc-Pk(nEl);  
 eq2(nEl) = q_lye0(nEl)*par.Const.CpLye*(T_in-T_k0(nEl)) + par.EL(nEl).nc*(U(nEl)-par.EL(nEl).Utn)*I(nEl) - ...
     (1/par.TherMo(nEl).Rt)*(T_k0(nEl)-par.EL(nEl).Ta);
 eq3(nEl) = U(nEl) - (par.U(nEl).r1+par.U(nEl).r2*T_k0(nEl))*I(nEl)/par.EL(nEl).A - par.U(nEl).s*log10(((par.U(nEl).t1+par.U(nEl).t2/T_k0(nEl)+...
     par.U(nEl).t3/(T_k0(nEl)^2))*I(nEl)/par.EL(nEl).A)+1) - par.EL(nEl).Urev; 
 end
 for nEl=1:par.N-1
     eq4(nEl) = U(nEl)*par.EL(nEl).nc - U(nEl+1)*par.EL(nEl+1).nc;
 end
 eq5 = (sum(Pk) - 63000);
 
 
%% Calculation of initial state vector
Feff = zeros(1,par.N);
nH2 = zeros(1,par.N);
qH2O_loss = zeros(1,par.N);
for nEl = 1:par.N
    Feff(nEl)= ((0.1*I(nEl)/par.EL(nEl).A)^2)/(par.U(nEl).f1+((0.1*I(nEl)/par.EL(nEl).A)^2))*par.U(nEl).f2;
    nH2(nEl) = Feff(nEl)*par.EL(nEl).nc*I(nEl)/(par.Const.ze*par.Const.FC);
    qH2O_loss(nEl) = nH2(nEl)*par.Const.Mwt;
end
flowH2 = nH2*0.0224136*3600
u_ini = U;
i_ini = I;
Tk_ini = T_k0;
for nEl = 1:par.N
    V_ini(nEl) = U(nEl)*(par.EL(nEl).nc);
end

Feff_ini = Feff;
nH2el_ini = nH2;
nH2out_ini = sum(nH2);
nO2 = nH2/2;
nO2el_ini = nO2;
nO2out_ini = sum(nO2);
%T_out = (q_lye0*T' - qH2O_loss*T')/(sum(q_lye0) - sum(qH2O_loss));  %temp after mixing of liquid streams from all electrolyzers
%Qcool = sum(q_lye0)*par.Const.CpLye*(T_out - T_in);                              %cooler duty at steady state, [J/s]
%SpecEl_ini = (Pk*(10^-3)*3600)./nH2el_ini;
% for i=1:N
% Qgen(i) = nc(i)*(U(i)-Utn(i))*I(i)
% end
%z0 = [u_ini, i_ini, Pk, Feff_ini, nH2el_ini, qH2O_loss, nH2out_ini, nO2el_ini, nO2out_ini, T_out];
%x0 = [Tk_ini, Qcool, V_ini, SpecEl_ini];

end
