clc
clear 
close all

%% Load parameters
N=3;
par = parElectrolyzer(N);

T = [80];
I = [0:0.1:5720]; %IdenMax = 220 mA/cm2

U=[];
Feff = [];

for nEl = 1:par.N  
    for k = 1:length(T)      
        %relation for Urev with T from LeRoy eqn. 58
        Urev(k) = 1.5184 - 1.5421e-3*(273+T(k)) + 9.523e-5*(273+T(k))*log((273+T(k))) + ...
            9.84e-8*(273+T(k))^2;
        
        u = Urev(k) + ((par.U(nEl).r1+par.U(nEl).r2*T(k)).*I./par.EL(nEl).A)...
            + (par.U(nEl).s*log10(((par.U(nEl).t1+par.U(nEl).t2/T(k)+...
            par.U(nEl).t3/(T(k)^2)).*I./par.EL(nEl).A)+1));
        
        feff = ((.1.*I./par.EL(nEl).A).^2)./(par.U(nEl).f1+((.1.*I./par.EL(nEl).A).^2))*par.U(nEl).f2;
        
        U=[U;u];
        Feff = [Feff;feff];     
    end   
end

rateH2 = Feff.*par.EL(nEl).nc.*I./(par.Const.ze*par.Const.FC); %in mol/s
vstd = 0.0224136;
nH2 = rateH2*vstd*3600; %in Nm3/hr

%% plot for UI characteristics

set(0, 'DefaultAxesFontSize', 12,  'DefaultLineLineWidth', 1.5)
fig1 = figure('NumberTitle', 'off', 'Name', 'UI characteristics');
hold on
% I-U curve
Iden = 0.1*I/par.EL(nEl).A;

plot(Iden,U(1,:),'b')
hold on
plot(Iden,U(2,:),'k')
hold on
plot(Iden,U(3,:),'r')

xlim([0, 220])
xlabel('Current Density, mA/cm^2')
ylabel('Voltage, V/cell')
grid on
legend('El 1 @ T = 80^{\circ}C','El 2 @ T = 80^{\circ}C','El 3 @ T = 80^{\circ}C')
legend('location','Northwest')
saveas(fig1,'UI_DegEl','epsc');

% %plot for U-I characteristics with the change in the temperature
% set(0, 'DefaultAxesFontSize', 12,  'DefaultLineLineWidth', 1.5)
% fig2 = figure('NumberTitle', 'off', 'Name', 'UI characteristics with temeperature');
% hold on
% % I-U curve
% Iden = 0.1*I/par.EL(nEl).A;
% 
% plot(Iden,U(1,:),'b')
% hold on
% plot(Iden,U(2,:),'k')
% hold on
% yline(1.229,'--b','LineWidth',1.5);
% yline(1.184,'--k','LineWidth',1.5);
% yline(1.482,'r','LineWidth',1.5);
% xlim([0, 220])
% xlabel('Current Density, mA/cm^2')
% ylabel('Voltage, V/cell')
% grid on
% legend('U @ T = 25^{\circ}C','U @ T = 80^{\circ}C','U_{rev} @ T = 25^{\circ}C',...
%     'U_{rev} @ T = 80^{\circ}C','U_{tn}')
% legend('location','Northwest')
% saveas(fig2,'UIChar','epsc')