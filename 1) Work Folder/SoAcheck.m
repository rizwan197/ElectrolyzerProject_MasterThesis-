clc
clear 
close all

T = 80;
par.Ncell = 230;
par.A = 2.6; %m2

Iden = 5:0.1:250;

par.F = 96485;
par.z = 2;
par.MwH2 = 2;
par.delG = 237*10^3; %J/mol

I = Iden.*par.A*10;

par.Urev = 1.5184 - 1.5421e-3*(273+T) + 9.523e-5*(273+T)*log((273+T)) + ...
    9.84e-8*(273+T)^2; %relation for Urev with T from LeRoy eqn. 58

par.U.r1 = 0.000218155;%vidar
% par.U.r1 = 3.799*10^-4; %using IPOPT
%par.U.r1 = 8.05*10^-5;%Ulleberg

par.U.r2 = -0.000000425; %vidar
% par.U.r2 = -2.5*10^-7;%Ulleberg

par.U.s = 0.1179375;%vidar
% par.U.s = 0.185;%ulleberg

par.U.t1 = -0.14529;%vidar
% par.U.t1 = -0.1443; %using IPOPT
%par.U.t1 = -0.1002;%ulleberg

par.U.t2 = 11.794;%vidar
% par.U.t2 = 8.424;%ulleberg

par.U.t3 = 395.68;%vidar
% par.U.t3 = 247.3;%ulleberg

par.U.f1 = 120;%vidar 
% par.U.f1 = 250;%ulleberg

par.U.f2 = 0.98;


u= par.Urev + (((par.U.r1 + par.U.r2*T).*I)./par.A) + par.U.s*log10(((par.U.t1+(par.U.t2/T)+...
    (par.U.t3/T^2)).*I/par.A)+1);
feff = (((0.1*I./par.A).^2)./(par.U.f1+(0.1*I./par.A).^2)).*par.U.f2;

nH2 = feff*par.Ncell.*I/(par.z*par.F);
VH2 = nH2*0.0224136*3600; %Nm3/hr

SpecEl = u.*I*par.Ncell./(1000.*VH2);
P = u.*I*par.Ncell/1000;

plot(Iden,u)
xlabel('Current density, mA/cm^2')
ylabel('Cell voltage, V/cell')
grid on

