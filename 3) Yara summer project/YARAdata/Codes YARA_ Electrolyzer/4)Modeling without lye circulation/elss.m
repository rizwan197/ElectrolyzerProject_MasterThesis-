function [x0,z0] = elss(param,n)
%This function file solves the steady state solution for each of the
%individual electrolyzer
%Here we will assume steady state temp of lye as 80C

%n = sequence of the electrolyzer

%Here we are solving following eqns for each electrolyzer:
%1)UI*nc-Power = 0;
%2)U - (((r1+r2*T)*I)/A) + s*log10(((t1+(t2/T)+(t3/T^2))*I/A)+1)-Urev = 0;

% Load CasADi
import casadi.*
global N

%% Load parameters
par = parElectrolyzer(N);

%parameter assignmemt
r1 = par.U(n).r1;
r2 = par.U(n).r2;
s = par.U(n).s;
t1 = par.U(n).t1;
t2 = par.U(n).t2;
t3 = par.U(n).t3;
f1 = par.U(n).f1;
f2 = par.U(n).f2;
         
Rt = par.TherMo(n).Rt;          

Utn = par.EL(n).Utn;            
Urev = par.EL(n).Urev;          
nc = par.EL(n).nc;              
A = par.EL(n).A;                
Ta = par.EL(n).Ta;              

ze = par.Const.ze;                   %number of electrons transferred per reaction 
FC = par.Const.FC;                   %faraday constant 


%% Define symbolic variables
x = SX.sym('x',2);
u=x(1);
i=x(2);

% Initial conditions for steady state calculations
Ps = param(n).Pel;                     %initial guess for power, [watts]
u0 = param(n).u0;                      %initial guess for voltage
i0 = param(n).i0;                      %initial guess for current
T = param(n).T0;                       %steady state temperature
% Ps = 4000;                            %initial guess for power, [watts]
% u0 = 1.3;                             %initial guess for voltage
% i0 = Ps/(nc*u0);                      %initial guess for current
                        
%% Solving steady state problem
g0 = u*i*nc-Ps;                                                         % P=UI*nc, i.e. Power to the electrolyzer, [Watts]
g1 = u - (r1+r2*T)*i/A - s*log10(((t1+t2/T+t3/(T^2))*i/A)+1) - Urev;    %U-I relationship

g = Function('g',{x},{[g0;g1]});
G = rootfinder('G','newton',g);
res = G([u0;i0]);

U = res(1);
I = res(2);

% Residuals check
% eq1 = U*I*nc-Ps
% eq2 = U - (r1+r2*T)*I/A - s*log10(((t1+t2/T+t3/(T^2))*I/A)+1) - Urev

%% Calculation of initial state vector
Feff=((0.1*I/A)^2)/(f1+((0.1*I/A)^2))*f2;
nH2=Feff*nc*I/(ze*FC);

u_ini = U;
i_ini = I;
Feff_ini = Feff;
nH2el_ini = nH2;
nH2out_ini = nH2;
nO2 = nH2/2;
nO2el_ini = nO2;
nO2out_ini = nO2;

%Calculation of Qgen and Qloss at steady state assuming steady state temp to be 80C
Qgen = nc*(U-Utn)*I;        %heat generated in the jth electrolyzer
Qloss = (T-Ta)/Rt;          %heat loss in the jth electrolyzer

z0 = [u_ini, i_ini, Feff_ini, nH2el_ini, nH2out_ini, nO2el_ini, nO2out_ini];
x0 = [Qgen, Qloss];

end
