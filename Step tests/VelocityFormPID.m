clc
clear
close all

% This is file cstrTpi.m
% INPUT: Time t and state vector y
% OUTPUT: derivatives DYDT
% Usage with ode15s:    [T,Y]=  ode15s('cstrTpi',[0 50],[2274 7726 444.0 0])

%% Load CasADi
addpath('/Users/mdrizwan/Documents/MATLAB/casadi-osx-matlabR2015a-v3.4.5')
import casadi.*

%% Define symbolic variables
x = SX.sym('x',3); 
p = SX.sym('p');
cA=x(1);
cB=x(2);
T=x(3);
xk = [2274;7726;444];%initial value

timeEnd = 20;       %final time
N = 1000;
dt = timeEnd/N;             %step size
ts = 0:dt:timeEnd;

%% Parameters
% All in SI units except time which is in minutes.
% I. Data (parameters and independent variables)
cAf=10000;             % feed concentration of A [mol/m3]
cBf=0;                 % feed concentration of B [mol/m3]
Tf=300;                % feed temperature [K]
V = 0.9;               % reactor volume [m3]
q = 0.12;              % volumetric flowrate [m3/min] (KAN ENDRES)
dhr =  -80e3;           % Heat of reaction [J/mol]
cp = 4.e3;              % heat capacity [J/kg]
rho = 1000;            % mean density [kg/m3]
U = 1000*60;           % heat transfer coefficient [J/min,m2,K]
A = 5;                 % heat transfer area [m2]
k1_400 = 0.1;          % rx1: rate constant at 400K [1/min]
k2_400 = 0.001;        % rx2: rate constant at 400K [1/min]
E1 = 60e3;             % rx1: activation energy [J/mol]
E2 = E1 - dhr;         % rx2: activation energy
R = 8.31;              % [J/K mol]
%% MV
u0 = 430;              % cooling temperature [K]
Tc = p;
%% PI controller
tauC    = 3;
k       = 0.5;
tau1    = 7;
Kc      = (1/k)*(tau1/tauC);
Ki      = Kc/tau1;

%% Model equations
% IIIa. Intermediate calculations
k1 = k1_400 * exp(-(E1/R) * (1/T - 1/400));
k2 = k2_400 * exp(-(E2/R) * (1/T - 1/400));
r = k1*cA - k2*cB;
Q = U*A*(Tc-T);

% IIIb. Evaluate derivatives of states
Vdcadt = q*cAf - q*cA - r*V;                   % [mol A/min]
Vdcbdt = q*cBf - q*cB + r*V;                   % [mol B/min]
mcpdTdt = rho*q*cp*(Tf-T) + r*V*(-dhr) + Q;    % [J/min]

f1 = Vdcadt/V;
f2 = Vdcbdt/V;
f3 = mcpdTdt / (rho*cp*V);

%% Model output and integrator
xdot=[f1; f2; f3];

%formulate the dynamics
ode = struct('x',x,'p',Tc,'ode',xdot);

%building integrator
opts = struct('tf',dt);         %integrating from 0 - dt 
F = integrator('F', 'cvodes', ode, opts);

%% Simulate the plant
err0 = 0;
x_n = xk;
MV=u0;
for i = 1:(length(ts)-1)
    
    %velocity form of the PI controller
    %velocity form calculates the change in the controller output
    err = 444-xk(3);
    u = u0 + Kc*(err-err0) + Ki*err*dt;
    err0 = err;
    u0=u;
    MV = [MV u];%for plotting
    
    % Integrate till the end of the interval
    Fk = F('x0',xk,'p',u);
    xk = full(Fk.xf);
    x_n = [x_n; full(Fk.xf)];
    
end
%%
T = x_n(3:3:end);
plot(ts,T)
hold on
plot(ts,MV)
ylim([430,445])

