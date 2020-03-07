function [x0,z0] = elss(param)
%This script solves for the steady state solution of the electrolyzer

% Load CasADi
import casadi.*

%% Load parameters
parElectrolyzer

%Here we are solving following steady state eqns:
%1)UI-34000 = 0;
%2)U - (((r1+r2*T)*I)/A) + s*log10(((t1+(t2/T)+(t3/T^2))*I/A)+1)-Urev = 0;
%3)aT-b = 0
%4)nH2 - (Feff*nc*I)/(ze*FC) = 0
%At steady state, nH2el = nH2out, nO2el = nO2out

%% Define symbolic variables
x = SX.sym('x',4);
u=x(1);
i=x(2);
T=x(3);
nH2=x(4);

%% Initial conditions for steady state calculations
Ps = param.Pel;                     %initial guess for power, [watts]
u0 = param.u0;                      %initial guess for voltage
i0 = param.i0;                      %initial guess for current
T0 = param.T0;                      %initial guess for electrolyser temperature, [celcius]
nH2o = param.nH2o;                  %initial guess for hydrogen flow rate, [mol/s]

%Ps = 3400;                         %initial guess for power
%u0 = 1.2;                          %initial guess for voltage
%i0 = Ps/u0;                        %initial guess for current
%T0 = 80;                           %initial electrolyser temperature, [celcius]

%% Solving steady state problem
g0 = u*i-Ps;% P=UI, i.e. Power to the electrolyzer, [Watts]
g1 = u - (r1+r2*T)*i/A - s*log10(((t1+t2/T+t3/T^2)*i/A)+1) - Urev;
g2 = ((1/tauT) + (Ccw/Ct)*(1-exp(-(hcond+hconv*i)/Ccw)))*T - ...
    (nc*(u-Utn)*1/Ct + (Ta/tauT) + (Ccw*Tcwi/Ct)*(1-exp(-(hcond+hconv*i)/Ccw)));
g3 = nH2 - (((.1*i/A)^2)/(f1+((.1*i/A)^2))*f2)*nc*i/(ze*FC);

g = Function('g',{x},{[g0;g1;g2;g3]});
G = rootfinder('G','newton',g);
res = G([u0;i0;T0;nH2o]);
U = res(1);
I = res(2);
T = res(3);
nH2 = res(4);

%% Calculation of initial state vector
Feff=((0.1*I/A)^2)/(f1+((0.1*I/A)^2))*f2;

u_ini = U;
i_ini = I;
a_ini = ((1/tauT) + (Ccw/Ct)*(1-exp(-(hcond+hconv*I)/Ccw)))*T;
b_ini = nc*(U-Utn)*1/Ct + (Ta/tauT) + (Ccw*Tcwi/Ct)*(1-exp(-(hcond+hconv*I)/Ccw));
Feff_ini = Feff;
nH2el_ini = nH2;
nH2out_ini = nH2;
nO2 = nH2/2;
nO2el_ini = nO2;
nO2out_ini = nO2;
Tel_ini = T;

z0 = [u_ini, i_ini, a_ini, b_ini, Feff_ini, nH2el_ini, nH2out_ini, nO2el_ini, nO2out_ini];
x0 = Tel_ini;
% kvlv=nH2/Vdisp;
end
