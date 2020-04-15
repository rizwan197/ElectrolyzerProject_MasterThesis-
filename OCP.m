function u1 = OCP(MPCinit,Tp,Tc,N)
import casadi.*
par = parElectrolyzer(N);
[lbx,lbz,lbu,ubx,ubz,ubu] = decision_var_bounds(par.N);

X0 = MPCinit.X0;     %X0 of the plant at t=0
dx0 = X0(6*par.N+6:end);
z0 = X0(1:6*par.N+5);
Xset = MPCinit.Xset;            %setpoint from the RTO
U0 = MPCinit.uk_prev;              %uk at t=0
h = Tp/Tc;
%% Define the Model equations
[x_var, z_var, p_var, alg, diff,~] = modelnew(N);
xk = [z_var;x_var];

%objective term
nu = 2*par.N+4;%length(p_var);
nz = 6*par.N+5;
nx = par.N+6;

Q = eye(nx+nz);
R = 0.5*eye(nu);
% L = (xk - Xset)'*Q*(xk-Xset) + (p_var-Uk_prev)'*R*(p_var-Uk_prev);
L = (xk - Xset)'*Q*(xk-Xset);
% Continuous time dynamics
f = Function('f',{x_var,z_var,p_var},{diff,alg,L},{'x','z','p'},{'xdot','zeval','qj'});

%% Direct Collocation
% Standard piece of code.

% Degree of interpolating polynomial
d = 3;

% Get collocation points
tau_root = [0, collocation_points(d, 'radau')];

% Coefficients of the collocation equation (of form xdot = C*x)
C = zeros(d+1,d+1);

% Coefficients of the continuity equation (Interpolating collocation points to get value at end of interval length [0-1])
D = zeros(d+1, 1);

% Coefficients of the quadrature function
B = zeros(d+1, 1);

% Construct polynomial basis
for j=1:d+1
    % Construct Lagrange polynomials to get the polynomial basis at the collocation point
    coeff = 1;
    for r=1:d+1
        if r ~= j
            coeff = conv(coeff, [1, -tau_root(r)]);
            coeff = coeff / (tau_root(j)-tau_root(r));
        end
    end
    % Evaluate the polynomial at the final time to get the coefficients of the continuity equation
    D(j) = polyval(coeff, 1.0);
    
    % Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the continuity equation
    pder = polyder(coeff);
    for r=1:d+1
        C(j,r) = polyval(pder, tau_root(r));
    end
    
    % Evaluate the integral of the polynomial to get the coefficients of the quadrature function
    pint = polyint(coeff);
    B(j) = polyval(pint, 1.0);
end

%% Build NLP solver
% empty nlp
w = {};
w0 = [];
lbw = [];
ubw = [];
J = 0;

g = {};
lbg = [];
ubg = [];

% initial conditions
X0 = MX.sym('X0',nx);
w = {w{:}, X0};
lbw = [lbw; dx0];
ubw = [ubw; dx0];
w0 = [w0; dx0];

% X0_par = MX.sym('X0_par',nx);
% % initial conditions
% g = {g{:},X0 - X0_par};
% lbg = [lbg;zeros(nx,1)];
% ubg = [ubg;zeros(nx,1)];

%% Building the NLP structure
Xk = X0;
% Uk_prev = U0;

for k = 0:Tp-1
    
    % Control Input (Uk)
    Uk = MX.sym(['U_' num2str(k)],nu);
    w = {w{:},Uk};
    lbw = [lbw;lbu];
    ubw = [ubw;ubu];
    w0 = [w0;U0];
    
    % Declaring New Collocation variables that are handles for solver
    Xkj = {};
    Zkj = {};
    
    for j = 1:d
        Xkj{j} = MX.sym(['X_' num2str(k) '_' num2str(j)],nx);
        Zkj{j} = MX.sym(['Z_' num2str(k) '_' num2str(j)],nz);
        
        w = {w{:},Xkj{j},Zkj{j}};
        lbw = [lbw;lbx;lbz];
        ubw = [ubw;ubx;ubz];
        w0 = [w0; dx0;z0];
    end
    %% Loop over collocation points
    Xk_end = D(1)*Xk;
    
    for j = 1:d
        % Expression for the state derivative of the collocation point j (collocation equation RHS i.e xdot = C*x)
        xp = C(1,j+1)*Xk;  % helper state
        for r = 1:d
            xp = xp + C(r+1,j+1)*Xkj{r};
        end
        
        %Calculating the diff,alg and quadrature @ collocation point
        [fj,zj,qj] =  f(Xkj{j},Zkj{j},Uk);
        
        % dynamic and algebraic constraint must be equal to zero
        g = {g{:},h*fj-xp,zj};
        lbg = [lbg;zeros(nx,1);zeros(nz,1)];
        ubg = [ubg;zeros(nx,1);zeros(nz,1)];
        
        % Add contribution to the end states (radau last point is at 1. for legendre we will need to interpolate)
        Xk_end = Xk_end + D(j+1)*Xkj{j};
    end
    %% Reached last collocation point
    %Adding the contribution to quadrature at the last collocation point
    J = J + (B(j+1)*qj*h) ;
    
    
%     % declaring additional constraints
%     uElconst = [];
%     for nEl=1:par.N-1
%         uElconst = [uElconst;z_var(nEl)-z_var(nEl+1)];
%     end
%     
%     Iden = MX.zeros(par.N,1);
%     for nEl = 1:par.N
%         Iden(nEl) = (0.1*z_var(par.N+nEl))/par.EL(nEl).A; %current density in mA/cm2
%     end
%     IdenMin = 32;   %minimum current density, 32 mA/cm2
%     IdenMax = 198.5;%maximum current density, 198.5 mA/cm2
%     
%     deltaT1 = x_var(par.N+5)-par.Tw_in;
%     deltaT2 = x_var(par.N+4)-x_var(par.N+6);
%     
%     Ptot = MX.zeros(1,1);
%     for nEl = 1:par.N
%         Ptot = Ptot + z_var(2*par.N+nEl);
%     end
%     
%     
%     g = {g{:},uElconst,Iden,deltaT1,deltaT2,Ptot};
%     lbg = [lbg;zeros(par.N-1,1);IdenMin*ones(par.N,1);2e-3;2e-3;0];
%     ubg = [ubg;zeros(par.N-1,1);IdenMax*ones(par.N,1);inf;inf;Pnet];
    
    
    % New NLP variable for state at end of interval
    Xk = MX.sym(['X_' num2str(k+1)], nx);
    w = {w{:},Xk};
    lbw = [lbw;lbx];
    ubw = [ubw;ubx];
    w0 = [w0; dx0];
    
    % Shooting Gap constraint
    g = {g{:},Xk_end-Xk};
    lbg = [lbg;zeros(nx,1)];
    ubg = [ubg;zeros(nx,1)];
    
end
%% Solving NLP
% create and solve NLP solver
opts = struct('warn_initial_bounds',false, ...
    'print_time',false, ...
    'ipopt',struct('print_level',1) ...
    );
nlp = struct('x',vertcat(w{:}),'f',J,'g',vertcat(g{:}));
solver = nlpsol('solver','ipopt',nlp,opts);

% Solve the NLP
sol = solver('x0', w0, 'lbx', lbw, 'ubx', ubw,...
    'lbg', lbg, 'ubg', ubg);
%% Collect the x,z,and u profiles
wMPC = full(sol.x);
for k = 1:nu
    uMPC(:,k) = wMPC(nx+k:(nu+d*(nx+nz)+nx):end);
end
u1 = uMPC(1,:);
for k= 1:nx
    xMPC(:,k) = wMPC(nx+nu+d*(nx+nz)+k:(nu+d*(nx+nz)+nx):end);
end
for k= 1:nz
    zMPC(:,k) = wMPC((nx+nu+d*nx+(d-1)*nz)+k:(nu+d*(nx+nz)+nx):end);
end
end