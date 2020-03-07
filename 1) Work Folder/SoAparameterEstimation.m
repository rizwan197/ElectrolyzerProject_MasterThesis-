
%This file solves the parameter values for r1 and t1

%Here we will assume steady state temp of lye as 80C

%n = sequence of the electrolyzer

%Here we are solving following eqns for each electrolyzer:
%1)SpecEl(supplier) = SpecEL(calculated) at Iden=220
%2)SpecEl(supplier) = SpecEL(calculated) at Iden=32

clc
clear
close all

% Load CasADi
import casadi.*


%% Load parameters
par.F = 96485;
par.z = 2;
par.Ncell = 230;
par.MwH2 = 2;
par.A = 2.6; %m2
par.delG = 237*10^3; %J/mol

T = 80;

par.Urev = 1.5184 - 1.5421e-3*(273+T) + 9.523e-5*(273+T)*log((273+T)) + ...
    9.84e-8*(273+T)^2; %relation for Urev with T from LeRoy eqn. 58


par.U.r2 = -2.5*10^-7;
par.U.s = 0.185;
par.U.t2 = 8.424;
par.U.t3 = 247.3;
par.U.f1 = 250;
par.U.f2 = 0.98;

%parameter assignmemt
%r1 = par.U(n).r1;
r2 = par.U.r2;
s = par.U.s;
%t1 = par.U(n).t1;
t2 = par.U.t2;
t3 = par.U.t3;
f1 = par.U.f1;
f2 = par.U.f2;
         
            
Urev = par.Urev;          
Ncell = par.Ncell;              
A = par.A;                
             
ze = par.z;                   %number of electrons transferred per reaction 
FC = par.F;                   %faraday constant 
vstd = 0.0224136; %m3/mol

%% Define symbolic variables
x = SX.sym('x',6);
r1=x(1);
t1=x(2);
umax=x(3);
umin=x(4);
Imax=x(5);
Imin=x(6);

% Initial conditions for steady state calculations
r10 = 8.05*10^-5;
t10 = -0.1002;
umax0 = 1.6639;
umin0 = 1.4305;
Imax0 = 199.2*A*10;
Imin0 = 35.56*A*10;

Pmax = 2134;
Pmin = 276.45;
VH2max = 485;
VH2min = 72.75;
                        
%% Solving steady state problem
feffmax = (((0.1*Imax./A).^2)./(f1+(0.1*Imax./A).^2)).*f2;
nH2max = feffmax*Ncell*Imax/(ze*FC);

feffmin = (((0.1*Imin./A).^2)./(f1+(0.1*Imin./A).^2)).*f2;
nH2min = feffmin*Ncell*Imin/(ze*FC);

g0 = Pmax - umax*Imax*Ncell/1000;
g1 = Pmin - umin*Imin*Ncell/1000;
g2 = umax - (Urev + (((r1 + r2*T).*Imax)./A) + s*log10(((t1+(t2/T)+(t3/T^2)).*Imax/A)+1));
g3 = umin - (Urev + (((r1 + r2*T).*Imin)./A) + s*log10(((t1+(t2/T)+(t3/T^2)).*Imin/A)+1));
g4 = VH2max - nH2max*0.0224136*3600; %Nm3/hr
g5 = VH2min - nH2min*0.0224136*3600; %Nm3/hr

g = Function('g',{x},{[g0;g1;g2;g3;g4;g5]});
G = rootfinder('G','newton',g);
res = G([r10;t10;umax0;umin0;Imax0;Imin0]);

r1 = full(res(1))
t1 = full(res(2))
umax = full(res(3))
umin = full(res(4))
Imax = full(res(5))
Imin = full(res(6))



