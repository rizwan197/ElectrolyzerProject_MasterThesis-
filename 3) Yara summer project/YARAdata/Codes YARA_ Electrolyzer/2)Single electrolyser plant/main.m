% Load CasADi
import casadi.*

%Using CasADi we are solving system of ODE and nonlinear algebric eqns simultaneously
%1)UI-Power = 0;
%2)U - (((r1+r2*T)*I)/A) - s*log10(((t1+(t2/T)+(t3/T^2))*I/A)+1) - Urev = 0;
%3)dT/dt = b-aT;
%4)a - ((1/tauT) + (Ccw/Ct)*(1-exp(-(hcond+hconv*I)/Ccw))) = 0;
%5)b - (nc*(U-Utn)*I/Ct + (Ta/tauT)+(Ccw*Tcwi/Ct)*(1-exp(-(hcond+hconv*I)/Ccw)))]) = 0
%6)Feff - ((.1*I/A)^2)/(f1+((.1*I/A)^2))*f2
%7)nH2 - Feff*nc*I/(ze*FC)
%8)nH2out - kvlvH2*VdispH2*sqrt(PstoH2-PoutH2)
%9)dPstoH2/dt = (TstoH2*Rg/VstoH2)*(nH2-nH2out)
%10)nO2out - kvlvO2*VdispO2*sqrt(PstoO2-PoutO2)
%11)dPstoO2/dt = (TstoO2*Rg/VstoO2)*(nO2-nO2out)

%% Loading parameters
parElectrolyzer

%% Calculation of valve constant
param.Pel = 2000;                   %power input to the electrolyzer, [Watts]
param.u0 = 1.2;                     %initial guess for voltage, [V]
param.i0 = param.Pel/param.u0;      %initial guess for current, [A]
param.T0 = 80;                      %initial guess for electrolyser temperature, [celcius]
param.nH2o = 0;                     %initial guess for hydrogen flow rate, [mol/s]

kvalveH2 = kvlv(param);             %valve constant for hydrogen outlet
kvalveO2 = kvlvO2(param);           %valve constant for oxygen outlet 

%% Steady state solution
[x0,z0] = elss(param);

nH2ss = z0(6);
nO2ss = z0(8);
T_ini = 51.7;                                               %temp at startup of the electrolyzer, [C] (manual intervention)
%T_ini = x0;
Psto_iniH2 = (nH2ss/(kvalveH2*VdispH2))^2 + PoutH2;         %initial H2 storage pressure (calculated from steady state solution)
%Psto_ini = Pout+40;
Psto_iniO2 = (nO2ss/(kvalveO2*VdispO2))^2 + PoutO2;

x0 = [T_ini Psto_iniH2 Psto_iniO2];

num_hr = 1;                                         %no. of hours
t0 = 1;                                             %start, [s)]
ts = 1;                                             %time step, [s]
tf = num_hr*60*60;                                  %final, [s]
tsamp = t0:ts:tf;
len = length(tsamp);                                %number of simulation time steps


%% Initialize variables
P = zeros(len,1);                       %power, [Watt]
P(1:100)=param.Pel;                     %incremental step change in power
P(101:end)=param.Pel/2;
P(1001:end)=param.Pel;
P(2001:end)=1.1*param.Pel;

Temp = zeros(len+1,1);                    %temp of the electrolyzer, [C]
Temp(1) = T_ini;

PstoH2 = zeros(len+1,1);                  %H2 storage pressure, [bar]
PstoH2(1) = full(Psto_iniH2);

PstoO2 = zeros(len+1,1);                  %O2 storage pressure, [bar]
PstoO2(1) = full(Psto_iniO2);

U = 1.2*ones(len,1);
I = zeros(len,1);
I_den = zeros(len,1);
Feff = zeros(len,1);
nH2in = zeros(len,1);
nH2out = zeros(len,1);
nO2in = zeros(len,1);
nO2out = zeros(len,1);

%% Solving for T_el,U,I,Feff,nH2,nO2
z = SX.sym('z',9); x = SX.sym('x',3); p = SX.sym('p',1);
dae = struct('x',x,'z',z,'p',p,'ode',[z(4)-z(3)*x(1);(TstoH2*Rg/VstoH2)*(z(6)-z(7));(TstoO2*Rg/VstoO2)*(z(8)-z(9))],...
    'alg',[z(1)*z(2)- p;...
    z(1) - (r1+r2*x(1))*z(2)/A - s*log10(((t1+t2/x(1)+t3/x(1)^2)*z(2)/A)+1) - Urev;...
    z(3) - ((1/tauT) + (Ccw/Ct)*(1-exp(-(hcond+hconv*z(2))/Ccw)));...
    z(4) - (nc*(z(1)-Utn)*z(2)/Ct + (Ta/tauT) + (Ccw*Tcwi/Ct)*(1-exp(-(hcond+hconv*z(2))/Ccw)));...
    z(5) - ((.1*z(2)/A)^2)/(f1+((.1*z(2)/A)^2))*f2;...
    z(6) - z(5)*nc*z(2)/(ze*FC);...
    z(7) - kvalveH2*VdispH2*sqrt(x(2)-PoutH2);...
    z(8) - z(6)/2;...
    z(9) - kvalveO2*VdispO2*sqrt(x(3)-PoutO2)]);
F = integrator('F', 'idas', dae);

for i=1:len
    %T=x(1);PstoH2=x(2);PstoO2=x(3);U=z(1);I=z(2);a=z(3);b=z(4);Feff=z(5);nH2=z(6);nH2out=z(7)
    r = F('x0',x0,'z0',z0,'p',P(i));
    x0 = full(r.xf);
    z0 = full(r.zf);
    
    %Calculation of compressor power
    PcompH2 = ((r.zf(6)*kvalveH2*R*Tel)/(alpha*(k-1)))*(((r.xf(2)/Pel)^((k-1)/k))-1);
    PcompO2 = ((r.zf(8)*kvalveO2*R*Tel)/(alpha*(k-1)))*(((r.xf(3)/Pel)^((k-1)/k))-1);%assuming same k and Tel for O2
    
    %Plotting variables
    Temp(i+1) = full(r.xf(1));
    PstoH2(i+1) = full(r.xf(2));
    PstoO2(i+1) = full(r.xf(3));
    U(i) = full(r.zf(1));
    I(i) = full(r.zf(2));
    I_den(i) = 0.1*I(i)/A;
    Feff(i) = full(r.zf(5));
    nH2in(i) = full(r.zf(6));
    nH2out(i) = full(r.zf(7));
    
    nO2in(i) = full(r.zf(8));
    nO2out(i) = full(r.zf(9));
    disp(i)
end

%% Plotting the results
figure(1)
subplot(3,1,1)
plot(I)
xlabel('Time, s')
ylabel('Current')
subplot(3,1,2)
plot(U)
xlabel('Time, s')
ylabel('Voltage')
subplot(3,1,3)
plot(I_den,U)
xlabel('Current Density, mA/cm^2')
ylabel('V/cell')

figure(2)
subplot(3,1,1)
plot(Temp)
xlabel('Time, s')
ylabel('Temperature')
subplot(3,1,2)
plot(PstoH2)
xlabel('Time, s')
ylabel('Storage pressure, H_2')
subplot(3,1,3)
plot(PstoO2)
xlabel('Time, s')
ylabel('Storage pressure, O_2')

figure(3)
plot(nH2in)
xlabel('Time, s')
ylabel('Flowrate, mol/s')
hold on
plot (nH2out)
legend('n_H_2 in','n_H_2 out')
xlabel('Time, s')
hold off

figure(4)
plot(nO2in)
xlabel('Time, s')
ylabel('Flowrate, mol/s')
hold on
plot (nO2out)
legend('n_O_2 in','n_O_2 out')
xlabel('Time, s')
hold off