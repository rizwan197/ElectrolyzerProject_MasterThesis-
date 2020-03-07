clc
clear all
close all

parU = struct('r1',8.05*10^-5,'r2',-2.5*10^-7,'s',0.185,'t1',-0.1002,'t2',8.424,'t3',247.3); %I-U curve parameters 
F = 96485;
z = 2;
Ncell = 21;
delG = 237*10^3; %kJ/mol
A = .25; %m2

T = [25 80];
for i=1:length(T)
    Urev(i) = 1.5184 - 1.5421e-3*(273+T(i)) + 9.523e-5*(273+T(i))*log((273+T(i))) + 9.84e-8*(273+T(i))^2; %relation for Urev with T from LeRoy eqn. 58
end
%Urev = [1.229 1.184]; 
I = [0:0.1:875];%A

%% I-U curve, Fig. 5 
U=[];
for i=1:length(T)
u = Urev(i) + (((parU.r1+parU.r2*T(i)).*I)./A) + parU.s*log10(((parU.t1+(parU.t2/T(i))+(parU.t3/T(i)^2)).*I/A)+1);
U = [U; u];  
end

set(0, 'DefaultAxesFontSize', 12,  'DefaultLineLineWidth',1.5)
fig1 = figure('NumberTitle', 'off', 'Name', 'UI characteristics');
hold on

plot(0.1*I/A,U(1,:),'b');% I-U curve of an electrolyser is in mA/cm2 vs V
hold on
plot(0.1*I/A,U(2,:),'k');
hold on

Urevplot = Urev'.*ones(1,length(0.1*I/A));
Utnplot = 1.482*ones(1,length(0.1*I/A));
plot(0.1*I/A,Urevplot(1,:),'b--')
hold on

plot(0.1*I/A,Urevplot(2,:),'k--')
hold on
plot(0.1*I/A,Utnplot,'r')
grid on

xlabel('Current Density, mA/cm^2')
ylabel('Voltage, V/cell')
hold on

legend('U @ T = 25', 'U @ T = 80','U_r_e_v @ T = 25','U_r_e_v @ T = 80', 'U_t_n @ T = 25, 80')
saveas(fig1,'UIchar','epsc');


% %nH2 vs current, Fig. 9
% parF = struct('f1',250,'f2',0.96);
% Feff = (((0.1*I./A).^2)./(parF.f1+(0.1*I./A).^2)).*parF.f2;
% rateH2 = Feff.*Ncell.*I./(z*F); %in mol/s
% vstd = 0.0224136;
% nH2 = rateH2*vstd*3600; 
% 
% figure(2)
% plot(I,nH2)
% xlabel('Current, A')
% ylabel('H_2 Flow Rate, Nm^3/hr')
% 
% %plot between efficiency percentage and current density, Fig. 6
% figure(3)
% plot(0.1*I/A,100*Feff)
% xlabel('Current Density, mA/cm^2')
% ylabel('Efficiency, %')
% 
% %plot between specific electricity and rate of H2 production, Vidar's plot
% %not current needs revision
% figure(4)
% %set(gca,'FontSize',14)
% E=(U.*I)*3600;%in Wh
% SpElec = 1e-6*E./rateH2;
% plot(rateH2,SpElec)
% ylim([40 80])
% ylabel('Specific electricity, MWh/mol H_2')
% xlabel('Molar flow rate of H_2, mol/hr')
% 
% figure(5)
% %set(gca,'FontSize',14)
% E=(U.*I)*3600;%in Wh
% SpElec = 1e-6*E./(rateH2);
% plot(rateH2,SpElec)
% ylim([40 80])
% ylabel('Specific electricity, MWh/mol H_2')
% xlabel('Molar flow rate of H_2, mol/hr')
% 
