clc
clear
close all

N = 3;
T = 80; 

par = parSoAEl(N,T); 
Imax = 2200*par.A;
I = [0:0.1:Imax]; %A, since max. current density = 220mA/cm^4

U=[];
Feff=[];

for i=1:N
    for k=1:length(T)
        u = par.Urev(k) + (((par.U(i).r1 + par.U(i).r2*T(k)).*I)./par.A) + par.U(i).s*log10(((par.U(i).t1+(par.U(i).t2/T(k))+...
            (par.U(i).t3/T(k)^2)).*I/par.A)+1);
        feff = (((0.1*I./par.A).^2)./(par.U(i).f1+(0.1*I./par.A).^2)).*par.U(i).f2;
        U=[U;u];
        Feff = [Feff;feff];
    end
end


nH2 = Feff.*par.Ncell.*I./(par.z*par.F); %in mol/s
vstd = 0.0224136;
rateH2 = nH2*vstd*3600; %in Nm3/hr

E=U.*I.*1e-6*par.Ncell; %in MW
kgH2 = nH2*par.MwH2*1e-6*3600; %in t/hr
specEl = (E./kgH2); %in MWh/tH2


nLHV = (rateH2*par.lhvH2)./(U.*I.*par.Ncell); %Electrolyzer efficiency
Es = par.lhvH2./1000*nLHV; %specific electricity consumption, kWh/Nm3

Pcomp = 1.5*(((rateH2*par.polyExp*par.R*par.Tel)./(3600*vstd*par.CompEff*(par.polyExp-1))).*...
    (((par.Psto/par.Pel)^((par.polyExp-1)/par.polyExp))-1))/1000; 
%Power consumed by compressor for storage of gases at steady state kWh/Nm3

Epl = Es + Pcomp;

%% plotting
set(0, 'DefaultAxesFontSize', 12,  'DefaultLineLineWidth', 1.5)
fig1 = figure('NumberTitle', 'off', 'Name', 'UI characteristics');
hold on
subplot(3,2,1)
plot((0.1*I./par.A),U(1,:),'b')
hold on
plot((0.1*I./par.A),U(2,:),'k')
hold on
plot((0.1*I./par.A),U(3,:),'r')
grid on
ylabel('Cell volatge, V/cell')
xlabel('Current Density, mA/cm^2')
legend('El 1','El 2','El 3')

subplot(3,2,2)
plot((0.1*I./par.A),nLHV(1,:),'b')
hold on
plot((0.1*I./par.A),nLHV(2,:),'k')
hold on
plot((0.1*I./par.A),nLHV(3,:),'r')
grid on
ylabel('Electrolyzer efficiency')
xlabel('Current Density, mA/cm^2')
legend('El 1','El 2','El 3')

subplot(3,2,3)
plot((0.1*I./par.A),Es(1,:),'b')
hold on
plot((0.1*I./par.A),Es(2,:),'k')
hold on
plot((0.1*I./par.A),Es(3,:),'r')
grid on
ylabel('Specific Electricity, kWh/Nm^3')
xlabel('Current Density, mA/cm^2')
legend('El 1','El 2','El 3')

subplot(3,2,4)
plot((0.1*I./par.A),E(1,:),'b')
hold on
plot((0.1*I./par.A),E(2,:),'k')
hold on
plot((0.1*I./par.A),E(3,:),'r')
grid on
ylabel('Power consumed, MW')
xlabel('Current Density, mA/cm^2')
legend('El 1','El 2','El 3')

subplot(3,2,5)
plot((0.1*I./par.A),rateH2(1,:),'b')
hold on
plot((0.1*I./par.A),rateH2(2,:),'k')
hold on
plot((0.1*I./par.A),rateH2(3,:),'r')
grid on
ylabel('H_2 produced, Nm^3 H_2/hr')
xlabel('Current Density, mA/cm^2')
legend('El 1','El 2','El 3')

% set(0, 'DefaultAxesFontSize', 12,  'DefaultLineLineWidth', 1.5)
% fig2 = figure('NumberTitle', 'off', 'Name', 'Balance of plant');
% hold on
subplot(3,2,6)
plot(rateH2(1,:),specEl(1,:),'b')
hold on
plot(rateH2(2,:),specEl(2,:),'k')
hold on
plot(rateH2(3,:),specEl(3,:),'r')
ylim([40 65])
grid on
ylabel('Specific Electricity (BOP), MWh/tH_2')
xlabel('H_2 produced, Nm^3 H_2/hr')
legend('El 1','El 2','El 3')