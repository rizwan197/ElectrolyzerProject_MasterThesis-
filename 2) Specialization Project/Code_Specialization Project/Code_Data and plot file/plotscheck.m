clc
clear 
close all

load data_Tinset7000_MVQcool_12hr

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
saveas(fig1,'TinsetTk_Ik_nH2k','epsc');

set(0, 'DefaultAxesFontSize', 10,  'DefaultLineLineWidth', 1.5)
fig2 = figure('NumberTitle', 'off', 'Name', 'Electrolyzer cell voltage');
hold on
plot(U)
xlabel('t, s')
ylabel('Cell voltage, V/cell')
legend('El 1','El 2', 'El 3')
grid on
saveas(fig2,'TinsetUk','epsc');

set(0, 'DefaultAxesFontSize', 10,  'DefaultLineLineWidth', 1.5)
fig3 = figure('NumberTitle', 'off', 'Name', 'Power Distribution');
hold on
subplot(2,4,5)
plot(P(:,1),'b')
xlabel('t, s')
ylabel('Power El_1, W')
grid on

subplot(2,4,6)
plot(P(:,2),'k')
xlabel('t, s')
ylabel('Power El_2, W')
grid on

subplot(2,4,7)
plot(P(:,3),'r')
xlabel('t, s')
ylabel('Power El_3, W')
grid on

subplot(2,4,8)
plot(P_net,'--')
xlabel('t, s')
ylabel('Total Power, W')
grid on

Qnet = Qlosslye+Qgen-Qloss;
Qtot=[];
for i=1:length(Qnet)
    Qtot(i) = sum(Qnet(i,:));
end
subplot(2,4,1)
plot(Qnet(:,1),'b')
xlabel('t, s')
ylabel('Qnet EL_1, Watts')
grid on
title('Electrolyzer 1')

subplot(2,4,2)
plot(Qnet(:,2),'k')
xlabel('t, s')
ylabel('Qnet EL_2, Watts')
grid on
title('Electrolyzer 2')

subplot(2,4,3)
plot(Qnet(:,3),'r')
xlabel('t, s')
ylabel('Qnet EL_3, Watts')
grid on
title('Electrolyzer 3')

subplot(2,4,4)
plot(Qtot)
xlabel('t, s')
ylabel('Qnet_t_o_t_a_l, Watts')
title('Overall assembly')
grid on
saveas(fig3,'TinsetQdisspNet','epsc');

set(0, 'DefaultAxesFontSize', 10,  'DefaultLineLineWidth', 1.5)
fig4 = figure('NumberTitle', 'off', 'Name', 'Energy dissipation in the electrolyzer assembly');
hold on
subplot(3,3,1)
plot(Qlosslye(:,1),'b')
xlabel('t, s')
ylabel('Qlye_o_u_t, Watts')
title('Electrolyzer 1')
grid on

subplot(3,3,2)
plot(Qlosslye(:,2),'k')
xlabel('t, s')
ylabel('Qlye_o_u_t, Watts')
title('Electrolyzer 2')
grid on

subplot(3,3,3)
plot(Qlosslye(:,3),'r')
xlabel('t, s')
ylabel('Qlye_o_u_t, Watts')
title('Electrolyzer 3')
grid on

subplot(3,3,4)
plot(Qgen(:,1),'b')
xlabel('t, s')
ylabel('Qgen, Watts')
grid on

subplot(3,3,5)
plot(Qgen(:,2),'k')
xlabel('t, s')
ylabel('Qgen, Watts')
grid on

subplot(3,3,6)
plot(Qgen(:,3),'r')
xlabel('t, s')
ylabel('Qgen, Watts')
grid on

subplot(3,3,7)
plot(Qloss(:,1),'b')
xlabel('t, s')
ylabel('Qloss_s_u_r_r, Watts')
grid on

subplot(3,3,8)
plot(Qloss(:,2),'k')
xlabel('t, s')
ylabel('Qloss_s_u_r_r, Watts')
grid on

subplot(3,3,9)
plot(Qloss(:,3),'r')
xlabel('t, s')
ylabel('Qloss_s_u_r_r, Watts')
grid on
saveas(fig4,'TinsetQdissp','epsc');

set(0, 'DefaultAxesFontSize', 10,  'DefaultLineLineWidth', 1.5)
fig5 = figure('NumberTitle', 'off', 'Name', 'Temperature of lye in recirculation loop');
hold on
subplot(2,1,2)
plot(Tout)
xlabel('t, s')
ylabel('T_o_u_t, ^{\circ}C')
grid on

subplot(2,1,1)
plot(Temp(2:end,1),'b')
hold on
plot(Temp(2:end,2),'k')
hold on
plot(Temp(2:end,3),'r')
xlabel('t, s')
ylabel('T_k, ^{\circ}C')
grid on
legend('El 1','El 2', 'El 3')

% subplot(2,2,1)
% plot(Qc)
% xlabel('t, s')
% ylabel('Q_c_o_o_l, watts')
% grid on
% 
% subplot(2,2,3)
% plot(PstoH2)
% xlabel('t, s')
% ylabel('Psto_H_2, bar')
% grid on
% 
% subplot(2,2,4)
% plot(PcompH2)
% xlabel('t, s')
% ylabel('P_c_o_m_p H2, Watts')
% grid on
saveas(fig5,'TinsetTout7000','epsc');

% set(0, 'DefaultAxesFontSize', 10,  'DefaultLineLineWidth', 1.5)
% fig6 = figure('NumberTitle', 'off', 'Name', 'Pressure in the gas storage tanks');
% hold on
% subplot(1,2,1)
% plot(PstoH2)
% xlabel('t, s')
% ylabel('Psto_H_2, bar')
% grid on
% 
% subplot(1,2,2)
% plot(PcompH2)
% xlabel('t, s')
% ylabel('P_c_o_m_p H2, Watts')
% grid on
% saveas(fig6,'zH2Psto','epsc');

%% plot of the step in MV

set(0, 'DefaultAxesFontSize', 10,  'DefaultLineLineWidth', 1.5)
fig7 = figure('NumberTitle', 'off', 'Name', 'Step change in the inlet lye flowrate');
hold on
% subplot(3,1,1)
% plot(qlye(:,1),'b')
% hold on
% plot(qlye(:,2),'k')
% hold on
% plot(qlye(:,3),'r')
% xlabel('t, s')
% ylabel('Q_l_y_e, g/s')
% grid on
% legend('El 1','El 2', 'El 3')

% subplot (3,1,2)
% plot(V_El,'r')
% xlabel('t, s')
% ylabel('V_E_l, V')
% grid on

subplot(2,1,1)
plot(Tin,'r')
hold on
plot(T_ini,'k--')
xlabel('t, s')
ylabel('T_i_n, ^{\circ}C')
legend('T_i_n','T_i_n_ _s_e_t')
grid on

subplot(2,1,2)
plot(Qc)
xlabel('t, s')
ylabel('Q_c_o_o_l, watts')
grid on

% subplot(4,1,1)
% plot(Tin,'r')
% xlabel('t, s')
% ylabel('z_H_2')
% grid on
saveas(fig7,'StepTinset7000','epsc');