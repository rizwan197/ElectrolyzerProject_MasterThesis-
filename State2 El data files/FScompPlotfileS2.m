clc
clear
close all

load('Data_CEl_S2_NewHex.mat')
load('Data_CEl_S2_OSHex.mat')
load('Data_CEl_S2_DegHex.mat')
load('Data_DCEl_S2_NewHex.mat')
load('Data_DCEl_S2_OSHex.mat')
load('Data_DCEl_S2_DegHex.mat')

load('Data_CEl_S2_NewHex_qlyefix')
load('Data_CEl_S2_OSHex_qlyefix')
load('Data_CEl_S2_DegHex_qlyefix')
load('Data_DCEl_S2_NewHex_qlyeFix.mat')
load('Data_DCEl_S2_OSHex_qlyeFix.mat')
load('Data_DCEl_S2_DegHex_qlyeFix.mat')


loss_F1S2_HXNew_qlvary = ((row_DC_S2_OS(:,end) - plot_C_S2_NewHex(:,end))./row_DC_S2_OS(:,end)).*100;       %Loss for F1S2_HXNew_qlvary w.r.t F2S2_HXOS_qlvary
loss_F1S2_HXNew_qlfix = ((row_DC_S2_OS(:,end) - plot_C_S2_NewHex_qlyefix(:,end))./row_DC_S2_OS(:,end)).*100;%Loss for F1S2_HXNew_qlfix w.r.t F2S2_HXOS_qlvary
loss_F1S2_HXOS_qlvary = ((row_DC_S2_OS(:,end) - plot_C_S2_OSHex(:,end))./row_DC_S2_OS(:,end)).*100;         %Loss for F1S2_HXOS_qlvary w.r.t F2S2_HXOS_qlvary
loss_F1S2_HXOS_qlfix = ((row_DC_S2_OS(:,end) - plot_C_S2_OSHex_qlyefix(:,end))./row_DC_S2_OS(:,end)).*100;  %Loss for F1S2_HXOS_qlfix w.r.t F2S2_HXOS_qlvary
loss_F1S2_HXDeg_qlvary = ((row_DC_S2_OS(:,end) - plot_C_S2_DegHex(:,end))./row_DC_S2_OS(:,end)).*100;       %Loss for F1S2_HXDeg_qlvary w.r.t F2S2_HXOS_qlvary
loss_F1S2_HXDeg_qlfix = ((row_DC_S2_OS(:,end) - plot_C_S2_DegHex_qlyefix(:,end))./row_DC_S2_OS(:,end)).*100;%Loss for F1S2_HXDeg_qlfix w.r.t F2S2_HXOS_qlvary

loss_F2S2_HXNew_qlvary = ((row_DC_S2_OS(:,end) - row_DC_S2_New(:,end))./row_DC_S2_OS(:,end)).*100;          %Loss for F2S2_HXNew_qlvary w.r.t F2S2_HXOS_qlvary
loss_F2S2_HXNew_qlfix = ((row_DC_S2_OS(:,end) - row_DC_S2_New_qlyeFix(:,end))./row_DC_S2_OS(:,end)).*100;   %Loss for F2S2_HXNew_qlfix w.r.t F2S2_HXOS_qlvary
loss_F2S2_HXOS_qlfix = ((row_DC_S2_OS(:,end) - row_DC_S2_OS_qlyeFix(:,end))./row_DC_S2_OS(:,end)).*100;     %Loss for F2S2_HXOS_qlfix w.r.t F2S2_HXOS_qlvary
loss_F2S2_HXDeg_qlvary = ((row_DC_S2_OS(:,end) - row_DC_S2_Deg(:,end))./row_DC_S2_OS(:,end)).*100;          %Loss for F2S2_HXDeg_qlvary w.r.t F2S2_HXOS_qlvary
loss_F2S2_HXDeg_qlfix = ((row_DC_S2_OS(:,end) - row_DC_S2_Deg_qlyeFix(:,end))./row_DC_S2_OS(:,end)).*100;   %Loss for F2S2_HXDeg_qlfix w.r.t F2S2_HXOS_qlvary


set(0, 'DefaultAxesFontSize', 12,  'DefaultLineLineWidth', 1.5)
fig1 = figure('NumberTitle', 'off');
hold on
plot(row_DC_S2_OS(:,1),loss_F1S2_HXNew_qlvary,'b');
hold on
plot(row_DC_S2_OS(:,1),loss_F1S2_HXOS_qlvary,'k');
hold on
plot(row_DC_S2_OS(:,1),loss_F1S2_HXDeg_qlvary,'r');
hold on
plot(row_DC_S2_OS(:,1),loss_F1S2_HXNew_qlfix,'--b');
hold on
plot(row_DC_S2_OS(:,1),loss_F1S2_HXOS_qlfix,'--k');
hold on
plot(row_DC_S2_OS(:,1),loss_F1S2_HXDeg_qlfix,'--r');

legend('S_2F_1HX_{New}q_{lye,var}','S_2F_1HX_{OS}q_{lye,var}','S_2F_1HX_{EoL}q_{lye,var}',...
    'S_2F_1HX_{New}q_{lye,fix}','S_2F_1HX_{OS}q_{lye,fix}','S_2F_1HX_{EoL}q_{lye,fix}')
legend('location','North')
xlabel('Input Power, [MW]')
ylabel('% Loss in production w.r.t S_2 F_2 HX_{OS} q_{lye,var}')
grid on
xlim([1,7])
ylim([-10,30])
saveas(fig1,'CompF1wrtS2F2opt','epsc');

set(0, 'DefaultAxesFontSize', 12,  'DefaultLineLineWidth', 1.5)
fig2 = figure('NumberTitle', 'off');
hold on
plot(row_DC_S2_OS(:,1),loss_F2S2_HXNew_qlvary,'b');
hold on
plot(row_DC_S2_OS(:,1),loss_F2S2_HXDeg_qlvary,'r');
hold on
plot(row_DC_S2_OS(:,1),loss_F2S2_HXNew_qlfix,'--b');
hold on
plot(row_DC_S2_OS(:,1),loss_F2S2_HXOS_qlfix,'--k');
hold on
plot(row_DC_S2_OS(:,1),loss_F2S2_HXDeg_qlfix,'--r');

legend('S_2F_2HX_{New}q_{lye,var}','S_2F_2HX_{EoL}q_{lye,var}',...
    'S_2F_2HX_{New}q_{lye,fix}','S_2F_2HX_{OS}q_{lye,fix}','S_2F_2HX_{EoL}q_{lye,fix}')
legend('location','North')
xlabel('Input Power, [MW]')
ylabel('% Loss in production w.r.t S_2 F_2 HX_{OS} q_{lye,var}')
grid on
xlim([1,7])
ylim([-10,30])
saveas(fig2,'CompF2wrtS2F2opt','epsc');