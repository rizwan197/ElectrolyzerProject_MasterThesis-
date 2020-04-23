clc
clear
close all

load('Data_CEl_S2_NewHex.mat')
load('Data_CEl_S2_OSHex.mat')
load('Data_CEl_S2_DegHex.mat')
load('Data_DCEl_S2_NewHex.mat')
load('Data_DCEl_S2_OSHex.mat')
load('Data_DCEl_S2_DegHex.mat')

load('Data_CEl_S2_DegHex_qlyefix')

loss_NewHex = ((row_DC_S2_New(:,end) - plot_C_S2_NewHex(:,end))./row_DC_S2_New(:,end)).*100;
loss_OSHex = ((row_DC_S2_OS(:,end) - plot_C_S2_OSHex(:,end))./row_DC_S2_OS(:,end)).*100;
loss_DegHex = ((row_DC_S2_Deg(:,end) - plot_C_S2_DegHex(:,end))./row_DC_S2_Deg(:,end)).*100;

set(0, 'DefaultAxesFontSize', 16,  'DefaultLineLineWidth', 1.5)
fig1 = figure('NumberTitle', 'off', 'Name', 'Loss between separate and shared BoP flowsheets for Hex design based on new electrolyzer');
hold on
plot(row_DC_S2_New(:,1),loss_NewHex,'b');
legend('q_l_y_e varying','q_l_y_e fixed')
xlabel('Input Power, [MW]')
ylabel('% Loss in production')
xlim([1,7])
ylim([-2,14])


set(0, 'DefaultAxesFontSize', 16,  'DefaultLineLineWidth', 1.5)
fig2 = figure('NumberTitle', 'off', 'Name', 'Loss between separate and shared BoP flowsheets for oversized Hex design');
hold on
plot(row_DC_S2_New(:,1),loss_OSHex,'k');
legend('q_l_y_e varying','q_l_y_e fixed')
xlabel('Input Power, [MW]')
ylabel('% Loss in production')
xlim([1,7])
ylim([-2,14])


set(0, 'DefaultAxesFontSize', 16,  'DefaultLineLineWidth', 1.5)
fig3 = figure('NumberTitle', 'off', 'Name', 'Loss between separate and shared BoP flowsheets for Hex design based on degraded electrolyzer');
hold on
plot(row_DC_S2_New(:,1),loss_DegHex,'r');
legend('q_l_y_e varying','q_l_y_e fixed')
xlabel('Input Power, [MW]')
ylabel('% Loss in production')
xlim([1,7])
ylim([-2,14])
