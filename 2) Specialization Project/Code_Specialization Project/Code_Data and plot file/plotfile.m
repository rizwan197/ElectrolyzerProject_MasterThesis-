clc
clear 
close all

load data_Vstep_3hr

set(0, 'DefaultAxesFontSize', 12,  'DefaultLineLineWidth', 1.5)
fig1 = figure('NumberTitle', 'off', 'Name', 'Electrolyzer current');
hold on
plot(I)
xlabel('Time, s')
ylabel('Current, A')
legend('El 1','El 2', 'El 3')
grid on
saveas(fig1,'VIk_3','epsc');


set(0, 'DefaultAxesFontSize', 12,  'DefaultLineLineWidth', 1.5)
fig2 = figure('NumberTitle', 'off', 'Name', 'Electrolyzer cell voltage');
hold on
plot(U)
xlabel('Time, s')
ylabel('Cell voltage, V/cell')
legend('El 1','El 2', 'El 3')
grid on
saveas(fig2,'VUk_3','epsc');

set(0, 'DefaultAxesFontSize', 12,  'DefaultLineLineWidth', 1.5)
fig3 = figure('NumberTitle', 'off', 'Name', 'Power Distribution');
hold on
subplot(4,1,1)
plot(P(:,1),'b')
xlabel('Time, s')
ylabel('Power El_1, W')
grid on

subplot(4,1,2)
plot(P(:,2),'k')
xlabel('Time, s')
ylabel('Power El_2, W')
grid on

subplot(4,1,3)
plot(P(:,3),'r')
xlabel('Time, s')
ylabel('Power El_3, W')
grid on

subplot(4,1,4)
plot(P_net)
xlabel('Time, s')
ylabel('Total Power, W')
grid on
saveas(fig3,'VPk_3','epsc');

set(0, 'DefaultAxesFontSize', 12,  'DefaultLineLineWidth', 1.5)
fig4 = figure('NumberTitle', 'off', 'Name', 'Hydrogen production rate');
hold on
plot(nH2elout)
hold on
plot(nH2in)
xlabel('Time, s')
ylabel('H_2 production rate, mol/s')
legend('El 1','El 2', 'El 3','H_2_n_e_t El')
grid on
saveas(fig4,'VnH2k_3','epsc');


set(0, 'DefaultAxesFontSize', 12,  'DefaultLineLineWidth', 1.5)
fig5 = figure('NumberTitle', 'off', 'Name', 'Net enthalpy change in the electrolyzer assembly');
hold on
Qnet = Qlosslye+Qgen-Qloss;
Qtot=[];
for i=1:length(Qnet)
    Qtot(i) = sum(Qnet(i,:));
end
plot(Qnet)
hold on
plot(Qtot)
xlabel('Time, s')
ylabel('Qnet, Watts')
legend('El_1','El_2','El_3','Q_t_o_t_a_l')
grid on
saveas(fig5,'VQdisspNet_3','epsc');

set(0, 'DefaultAxesFontSize', 12,  'DefaultLineLineWidth', 1.5)
fig6 = figure('NumberTitle', 'off', 'Name', 'Energy dissipation in the electrolyzer assembly');
hold on
subplot(3,3,1)
plot(Qlosslye(:,1))
xlabel('Time, s')
ylabel('Qlye_o_u_t, Watts')
title('Electrolyzer 1')
grid on

subplot(3,3,2)
plot(Qlosslye(:,2))
xlabel('Time, s')
ylabel('Qlye_o_u_t, Watts')
title('Electrolyzer 2')
grid on

subplot(3,3,3)
plot(Qlosslye(:,3))
xlabel('Time, s')
ylabel('Qlye_o_u_t, Watts')
title('Electrolyzer 3')
grid on

subplot(3,3,4)
plot(Qgen(:,1))
xlabel('Time, s')
ylabel('Qgen, Watts')
grid on

subplot(3,3,5)
plot(Qgen(:,2))
xlabel('Time, s')
ylabel('Qgen, Watts')
grid on

subplot(3,3,6)
plot(Qgen(:,3))
xlabel('Time, s')
ylabel('Qgen, Watts')
grid on

subplot(3,3,7)
plot(Qloss(:,1))
xlabel('Time, s')
ylabel('Qloss_s_u_r_r, Watts')
grid on

subplot(3,3,8)
plot(Qloss(:,2))
xlabel('Time, s')
ylabel('Qloss_s_u_r_r, Watts')
grid on

subplot(3,3,9)
plot(Qloss(:,3))
xlabel('Time, s')
ylabel('Qloss_s_u_r_r, Watts')
grid on

saveas(fig6,'VQdissp_3','epsc');

set(0, 'DefaultAxesFontSize', 12,  'DefaultLineLineWidth', 1.5)
fig7 = figure('NumberTitle', 'off', 'Name', 'Electrolyzer temperature');
hold on
plot(Temp(2:end,:))
xlabel('Time, s')
ylabel('T_k, ^{\circ}C')
legend('El 1','El 2', 'El 3')
grid on
saveas(fig7,'VTk_3','epsc');

set(0, 'DefaultAxesFontSize', 12,  'DefaultLineLineWidth', 1.5)
fig8 = figure('NumberTitle', 'off', 'Name', 'Temperature of lye in recirculation loop');
hold on
subplot(2,1,2)
plot(Tout)
xlabel('Time, s')
ylabel('T_o_u_t, ^{\circ}C')
grid on

subplot(2,1,1)
plot(Tin)
xlabel('Time, s')
ylabel('T_i_n, ^{\circ}C')
grid on
saveas(fig8,'VTout_3','epsc');

%% plot of the step in MV

set(0, 'DefaultAxesFontSize', 12,  'DefaultLineLineWidth', 1.5)
fig9 = figure('NumberTitle', 'off', 'Name', 'Step change in the inlet lye flowrate');
hold on
plot(qlye)
xlabel('Time, s')
ylabel('Q_l_y_e, g/s')
grid on
legend('El 1','El 2', 'El 3')
saveas(fig9,'StepV_3','epsc');

% set(0, 'DefaultAxesFontSize', 12,  'DefaultLineLineWidth', 1.5)
% fig9 = figure('NumberTitle', 'off', 'Name', 'Step change in the Electrolyzer voltage');
% hold on
% plot(V_El)
% xlabel('Time, s')
% ylabel('V_E_l, V')
% grid on
% saveas(fig9,'StepV_3','epsc');

% set(0, 'DefaultAxesFontSize', 12,  'DefaultLineLineWidth', 1.5)
% fig9 = figure('NumberTitle', 'off', 'Name', 'Step change in cooling duty');
% hold on
% plot(Qc)
% xlabel('Time, s')
% ylabel('Q_c_o_o_l, watts')
% grid on
% saveas(fig9,'StepQc_3','epsc');