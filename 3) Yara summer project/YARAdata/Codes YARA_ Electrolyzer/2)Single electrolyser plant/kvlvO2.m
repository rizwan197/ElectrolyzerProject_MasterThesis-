function kvlv = kvlvO2(param)
%This script calculates the valve constant for O2 outlet.

% Load CasADi
import casadi.*

%% Load parameters
parElectrolyzer

%Here we are solving following steady state eqns:
%1)UI-34000 = 0;
%2)U - (((r1+r2*T)*I)/A) + s*log10(((t1+(t2/T)+(t3/T^2))*I/A)+1)-Urev = 0;
%3)aT-b = 0

%% Define symbolic variables
x = SX.sym('x',3);
u=x(1);
i=x(2);
T=x(3);

%% Initial conditions for steady state calculations
Ps = param.Pel;                     %initial guess for power, [watts]
u0 = param.u0;                      %initial guess for voltage
i0 = param.i0;                      %initial guess for current
T0 = param.T0;                      %initial electrolyser temperature, [celcius]

%Ps = 3400;                         %initial guess for power
%u0 = 1.2;                          %initial guess for voltage
%i0 = Ps/u0;                        %initial guess for current
%T0 = 80;                           %initial electrolyser temperature, [celcius]

%% Solving steady state problem
g0 = u*i-Ps;% P=UI, i.e. Power to the electrolyzer in Watts
g1 = u - (r1+r2*T)*i/A - s*log10(((t1+t2/T+t3/T^2)*i/A)+1) - Urev;
g2 = ((1/tauT) + (Ccw/Ct)*(1-exp(-(hcond+hconv*i)/Ccw)))*T - ...
    (nc*(u-Utn)*1/Ct + (Ta/tauT) + (Ccw*Tcwi/Ct)*(1-exp(-(hcond+hconv*i)/Ccw)));

g = Function('g',{x},{[g0;g1;g2]});
G = rootfinder('G','newton',g);
res = G([u0;i0;T0]);
%U=res(1);
I=res(2);
%T=res(3);

Feff=((0.1*I/A)^2)/(f1+((0.1*I/A)^2))*f2;
nH2=(Feff*nc*I)/(ze*FC);
nO2=nH2/2;
kvlv=nO2/VdispO2;
end
