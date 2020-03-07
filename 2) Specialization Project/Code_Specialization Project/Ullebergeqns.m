clc
clear
close all

N=3;
T = [70 80];  
I = [0:0.1:875]; %A, since max. current density = 350mA/cm^4
par = parEL(N,T);

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


rateH2 = Feff.*par.Ncell.*I./(par.z*par.F); %in mol/s
vstd = 0.0224136;
nH2 = rateH2*vstd*3600; %in Nm3/hr

E=U.*I.*1e-6*par.Ncell; %in MW
kgH2 = rateH2*par.MwH2*1e-6*3600; %in t/hr
specEl = (E./kgH2); %in MWh/tH2


%% plot for UI characteristics

set(0, 'DefaultAxesFontSize', 12,  'DefaultLineLineWidth', 1.5)
fig1 = figure('NumberTitle', 'off', 'Name', 'UI characteristics');
hold on
% I-U curve
Iden = 0.1*I/par.A;

plot(Iden,U(1,:),'b')
hold on
plot(Iden,U(3,:),'k')
hold on
plot(Iden,U(5,:),'r')
hold on
plot(Iden,U(2,:),'b--')
hold on
plot(Iden,U(4,:),'k--')
hold on
plot(Iden,U(6,:),'r--')
grid on
xlabel('Current Density, mA/cm^2')
ylabel('Voltage, V/cell')
legend('El 1 @ T=70','El 2 @ T=70','El 3 @ T=70','El 1 @ T=80','El 2 @ T=80','El 3 @ T=80')
saveas(fig1,'UITemp','epsc');


%% plot between specific electricity and rate of H2 production
set(0, 'DefaultAxesFontSize', 12,  'DefaultLineLineWidth', 1.5)
fig2 = figure('NumberTitle', 'off', 'Name', 'Specific Electricity');
hold on
plot(kgH2(1,:),specEl(1,:),'b')
hold on
plot(kgH2(3,:),specEl(3,:),'k')
hold on
plot(kgH2(5,:),specEl(5,:),'r')
hold on
plot(kgH2(2,:),specEl(2,:),'b--')
hold on
plot(kgH2(4,:),specEl(4,:),'k--')
hold on
plot(kgH2(6,:),specEl(6,:),'r--')

ylim([40 65])
grid on

legend('El 1 @ T=70','El 2 @ T=70','El 3 @ T=70','El 1 @ T=80','El 2 @ T=80','El 3 @ T=80')
ylabel('Specific Electricity, MWh/tH_2')
xlabel('flowrate, tH_2/hr')
saveas(fig2,'Weff','epsc');