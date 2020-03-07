clc
clear

%% Load CasADi
import casadi.*

%% Load parameters
FC = 96485;
ze = 2;
Ncell = 230;
par.MwH2 = 2;
A = 2.6; %m2
par.delG = 237*10^3; %J/mol
T = 80;
Urev = 1.5184 - 1.5421e-3*(273+T) + 9.523e-5*(273+T)*log((273+T)) + ...
    9.84e-8*(273+T)^2; %relation for Urev with T from LeRoy eqn. 58
r2 = -2.5*10^-7;
s = 0.185;
t2 = 8.424;
t3 = 247.3;
f1 = 250;
f2 = 0.98;
vstd = 0.0224136; %m3/mol

%% Defining System
% declaring variables
x = SX.sym('x',6);

r1=x(1);
t1=x(2);
umax=x(3);
umin=x(4);
Imax=x(5);
Imin=x(6);

% declaring parameters
Pmax = 2134;
Pmin = 276.45;
VH2max = 485;
VH2min = 72.75;
r10 = 8.05*10^-5;
t10 = -0.1002;

% declaring constraints
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

% objective term
L = 10^5*(r1-r10)^2 + (t1-t10)^2;

% Optimization problem
% preparing symbolic variables
w = {};
% preparing numeric variables and bounds
w0 = [];
lbw = [];
ubw = [];

% declaring them symbolic
w = {w{:},r1,t1,umax,umin,Imax,Imin};
% declaring them numerically
w0 = [w0;8.05*10^-5;-0.1002;1.6221;1.4447;220;32]; 
%initial guess bounds
lbw = [lbw;-inf;-inf;-inf;-inf;-inf;-inf];
ubw = [ubw;inf;inf;inf;inf;inf;inf]; 

% preparing symbolic constraints
g = {};
% preparing numeric bounds
lbg = [];
ubg = [];

% declaring constraints
g = {g{:},g0,g1,g2,g3,g4,g5};
lbg = [lbg;0;0;0;0;0;0];
ubg = [ubg;0;0;0;0;0;0];

% optimization objective function
% By default, casadi always minimizes the problem. 
% Since we want to minimize it, we have to write:
J = L; 

% formalize it into an NLP problem
% the main difference is the fourth argument 'p'!
nlp = struct('x',vertcat(w{:}),'g',vertcat(g{:}),'f',J);

% assign solver - IPOPT in this case
solver = nlpsol('solver','ipopt',nlp);


% solve - using the previous defined initial guess and bounds
sol = solver('x0',w0,'lbx',lbw,'ubx',ubw,'lbg',lbg,'ubg',ubg);

% Extract Solution
% we use the function full() to transform the variables from Casadi type 
% to 'actual' numbers

% variables values at optimum
w_opt = full(sol.x) 
% objective function
costOF = full(sol.f); 
% constraints values
M_opt = full(sol.g); 
% multipliers
lam_M = full(sol.lam_g); 
% parameters multipliers
lam_theta = full(sol.lam_p); 