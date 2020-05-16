clc
clear
close all

load('Data_CEl_S1_NewHex.mat')
load('Data_CEl_S1_OSHex.mat')
load('Data_CEl_S1_DegHex.mat')
load('Data_DCEl_S1_NewHex.mat')
load('Data_DCEl_S1_OSHex.mat')
load('Data_DCEl_S1_DegHex.mat')

load('Data_CEl_S1_NewHex_qlyefix')
load('Data_CEl_S1_OSHex_qlyefix')
load('Data_CEl_S1_DegHex_qlyefix')
load('Data_DCEl_S1_NewHex_qlyeFix.mat')
load('Data_DCEl_S1_OSHex_qlyeFix.mat')
load('Data_DCEl_S1_DegHex_qlyeFix.mat')


loss_F1S1_HXNew_qlvary = ((row_DC_S1_OS(:,end) - plot_C_S1_NewHex(:,end))./row_DC_S1_OS(:,end)).*100;       %Loss for F1S1_HXNew_qlvary w.r.t F2S1_HXOS_qlvary
loss_F1S1_HXNew_qlfix = ((row_DC_S1_OS(:,end) - plot_C_S1_NewHex_qlyefix(:,end))./row_DC_S1_OS(:,end)).*100;%Loss for F1S1_HXNew_qlfix w.r.t F2S1_HXOS_qlvary
loss_F1S1_HXOS_qlvary = ((row_DC_S1_OS(:,end) - plot_C_S1_OSHex(:,end))./row_DC_S1_OS(:,end)).*100;         %Loss for F1S1_HXOS_qlvary w.r.t F2S1_HXOS_qlvary
loss_F1S1_HXOS_qlfix = ((row_DC_S1_OS(:,end) - plot_C_S1_OSHex_qlyefix(:,end))./row_DC_S1_OS(:,end)).*100;  %Loss for F1S1_HXOS_qlfix w.r.t F2S1_HXOS_qlvary
loss_F1S1_HXDeg_qlvary = ((row_DC_S1_OS(:,end) - plot_C_S1_DegHex(:,end))./row_DC_S1_OS(:,end)).*100;       %Loss for F1S1_HXDeg_qlvary w.r.t F2S1_HXOS_qlvary
loss_F1S1_HXDeg_qlfix = ((row_DC_S1_OS(:,end) - plot_C_S1_DegHex_qlyefix(:,end))./row_DC_S1_OS(:,end)).*100;%Loss for F1S1_HXDeg_qlfix w.r.t F2S1_HXOS_qlvary

loss_F2S1_HXNew_qlvary = ((row_DC_S1_OS(:,end) - row_DC_S1_New(:,end))./row_DC_S1_OS(:,end)).*100;          %Loss for F2S1_HXNew_qlvary w.r.t F2S1_HXOS_qlvary
loss_F2S1_HXNew_qlfix = ((row_DC_S1_OS(:,end) - row_DC_S1_New_qlyeFix(:,end))./row_DC_S1_OS(:,end)).*100;   %Loss for F2S1_HXNew_qlfix w.r.t F2S1_HXOS_qlvary
loss_F2S1_HXOS_qlfix = ((row_DC_S1_OS(:,end) - row_DC_S1_OS_qlyeFix(:,end))./row_DC_S1_OS(:,end)).*100;     %Loss for F2S1_HXOS_qlfix w.r.t F2S1_HXOS_qlvary
loss_F2S1_HXDeg_qlvary = ((row_DC_S1_OS(:,end) - row_DC_S1_Deg(:,end))./row_DC_S1_OS(:,end)).*100;          %Loss for F2S1_HXDeg_qlvary w.r.t F2S1_HXOS_qlvary
loss_F2S1_HXDeg_qlfix = ((row_DC_S1_OS(:,end) - row_DC_S1_Deg_qlyeFix(:,end))./row_DC_S1_OS(:,end)).*100;   %Loss for F2S1_HXDeg_qlfix w.r.t F2S1_HXOS_qlvary


figure()
plot(row_DC_S1_OS(:,1),loss_F1S1_HXNew_qlvary,'b');
hold on
plot(row_DC_S1_OS(:,1),loss_F1S1_HXOS_qlvary,'k');
hold on
plot(row_DC_S1_OS(:,1),loss_F1S1_HXDeg_qlvary,'r');
hold on
plot(row_DC_S1_OS(:,1),loss_F2S1_HXNew_qlvary,'--b');
hold on
plot(row_DC_S1_OS(:,1),loss_F2S1_HXDeg_qlvary,'--r');

legend('F_1S_1HX_{New}q_{lye,var}','F_1S_1HX_{OS}q_{lye,var}','F_1S_1HX_{Deg}q_{lye,var}','F_2S_1HX_{New}q_{lye,var}','F_2S_1HX_{Deg}q_{lye,var}')
xlabel('Input Power, [MW]')
ylabel('% Loss in production w.r.t F_2S_1HX_{OS}q_{lye,var}')
xlim([1,7])

figure()
plot(row_DC_S1_OS(:,1),loss_F1S1_HXNew_qlfix,'b');
hold on
plot(row_DC_S1_OS(:,1),loss_F1S1_HXOS_qlfix,'k');
hold on
plot(row_DC_S1_OS(:,1),loss_F1S1_HXDeg_qlfix,'r');
hold on
plot(row_DC_S1_OS(:,1),loss_F2S1_HXNew_qlfix,'--b');
hold on
plot(row_DC_S1_OS(:,1),loss_F2S1_HXOS_qlfix,'--k');
hold on
plot(row_DC_S1_OS(:,1),loss_F2S1_HXDeg_qlfix,'--r');

legend('F_1S_1HX_{New}q_{lye,fix}','F_1S_1HX_{OS}q_{lye,fix}','F_1S_1HX_{Deg}q_{lye,fix}','F_2S_1HX_{New}q_{lye,fix}','F_2S_1HX_{OS}q_{lye,fix}','F_2S_1HX_{Deg}q_{lye,fix}')
xlabel('Input Power, [MW]')
ylabel('% Loss in production w.r.t F_2S_1HX_{OS}q_{lye,var}')
xlim([1,7])

