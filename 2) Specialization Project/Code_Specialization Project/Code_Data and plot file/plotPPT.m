clc
clear 
close all

load data_qlye1step3000_MVQcool_12hr

set(0, 'DefaultAxesFontSize', 10,  'DefaultLineLineWidth', 1.5)
fig1 = figure('NumberTitle', 'off', 'Name', 'Electrolyzer temperature');
hold on
subplot(3,3,1)
plot(Temp(2:end,1),'b')
xlabel('t, s')
ylabel('T_k, ^{\circ}C')
title('Electrolyzer 1')
grid on

subplot(3,3,2)
plot(Temp(2:end,2),'k')
xlabel('t, s')
ylabel('T_k, ^{\circ}C')
title('Electrolyzer 2')
grid on

subplot(3,3,3)
plot(Temp(2:end,3),'r')
xlabel('t, s')
ylabel('T_k, ^{\circ}C')
title('Electrolyzer 3')
grid on

subplot(3,3,4)
plot(I(:,1),'b')
xlabel('t, s')
ylabel('Current, A')
grid on

subplot(3,3,5)
plot(I(:,2),'k')
xlabel('t, s')
ylabel('Current, A')
grid on

subplot(3,3,6)
plot(I(:,3),'r')
xlabel('t, s')
ylabel('Current, A')
grid on

subplot(3,3,7)
plot(nH2elout(:,1),'b')
xlabel('t, s')
ylabel('n H_2, mol/s')
grid on

subplot(3,3,8)
plot(nH2elout(:,2),'k')
xlabel('t, s')
ylabel('n H_2, mol/s')
grid on

subplot(3,3,9)
plot(nH2elout(:,3),'r')
xlabel('t, s')
ylabel('n H_2, mol/s')
grid on
saveas(fig1,'PPTqlye1SRTk_Ik_nH2k','epsc');

set(0, 'DefaultAxesFontSize', 10,  'DefaultLineLineWidth', 1.5)
fig5 = figure('NumberTitle', 'off', 'Name', 'Temperature of lye in recirculation loop');
hold on
subplot(2,2,2)
plot(Tout)
xlabel('t, s')
ylabel('T_o_u_t, ^{\circ}C')
grid on

subplot(2,2,1)
plot(Tin)
xlabel('t, s')
ylabel('T_i_n, ^{\circ}C')
ylim([64.5 66.5])
grid on

subplot(2,2,3)
plot(PstoH2)
xlabel('t, s')
ylabel('Psto_H_2 bar')
grid on

subplot(2,2,4)
plot(PcompH2)
xlabel('t, s')
ylabel('P_c_o_m_p_H_2, Watts')
grid on
saveas(fig5,'PPTqlye1SRAux','epsc');

%% plot of the step in MV

set(0, 'DefaultAxesFontSize', 10,  'DefaultLineLineWidth', 1.5)
fig7 = figure('NumberTitle', 'off', 'Name', 'Step change in the inlet lye flowrate');
hold on
subplot(2,1,1)
plot(qlye(:,1),'b')
hold on
plot(qlye(:,2),'k')
hold on
plot(qlye(:,3),'r')
xlabel('t, s')
ylabel('Q_l_y_e, g/s')
grid on
legend('El 1','El 2', 'El 3')

% subplot (4,1,3)
% plot(V_El,'r')
% xlabel('t, s')
% ylabel('V_E_l, V')
% grid on

subplot(2,1,2)
plot(Qc,'r')
xlabel('t, s')
ylabel('Q_c_o_o_l, watts')
grid on

% subplot(4,1,1)
% plot(ZH2,'r')
% xlabel('t, s')
% ylabel('z_H_2')
% grid on
% 
% subplot(5,1,5)
% plot(ZO2,'r')
% xlabel('t, s')
% ylabel('z_O_2')
% grid on

saveas(fig7,'PPTStepQlye1SR','epsc');