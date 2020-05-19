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

loss_NewHex = ((row_DC_S2_New(:,end) - plot_C_S2_NewHex(:,end))./row_DC_S2_New(:,end)).*100;
loss_OSHex = ((row_DC_S2_OS(:,end) - plot_C_S2_OSHex(:,end))./row_DC_S2_OS(:,end)).*100;
loss_DegHex = ((row_DC_S2_Deg(:,end) - plot_C_S2_DegHex(:,end))./row_DC_S2_Deg(:,end)).*100;

%plots section 3.1
set(0, 'DefaultAxesFontSize', 12,  'DefaultLineLineWidth', 1.5)
fig1 = figure('NumberTitle', 'off');
hold on
plot(row_DC_S2_New(:,1),loss_NewHex,'b');
hold on
plot(row_DC_S2_New(:,1),loss_OSHex,'k');
hold on
plot(row_DC_S2_New(:,1),loss_DegHex,'r');
legend('S_2F_1HX_{New}q_{lye,var} vs S_2F_2HX_{New}q_{lye,var}',...
    'S_2F_1HX_{OS}q_{lye,var} vs S_2F_2HX_{OS}q_{lye,var}',...
    'S_2F_1HX_{Deg}q_{lye,var} vs S_2F_2HX_{Deg}q_{lye,var}')
legend('location','North')
xlabel('Input Power, [MW]')
ylabel('% Loss in production')
grid on
xlim([1,7])
ylim([-2,14])
saveas(fig1,'FScomp','epsc');

%plots for section 3.2
%loss w.r.t OSHex design for CEl and DCEl with varying lye flowrate
loss_CEl_S2_NewHexComp = ((plot_C_S2_OSHex(:,end) - plot_C_S2_NewHex(:,end))./(plot_C_S2_OSHex(:,end))).*100;
loss_CEl_S2_DegHexComp = ((plot_C_S2_OSHex(:,end) - plot_C_S2_DegHex(:,end))./(plot_C_S2_OSHex(:,end))).*100;
loss_DCEl_S2_NewHexComp = ((row_DC_S2_OS(:,end) - row_DC_S2_New(:,end))./row_DC_S2_OS(:,end)).*100;
loss_DCEl_S2_DegHexComp = ((row_DC_S2_OS(:,end) - row_DC_S2_Deg(:,end))./row_DC_S2_OS(:,end)).*100;

% set(0, 'DefaultAxesFontSize', 1,  'DefaultLineLineWidth', 0.3)
% fig2 = figure('NumberTitle', 'off', 'Name', 'Loss of Hex design vs OS Hex design');
% hold on
figure()
subplot(2,1,1)
plot(row_DC_S2_New(:,1),loss_CEl_S2_NewHexComp,'b');
hold on
plot(row_DC_S2_New(:,1),loss_CEl_S2_DegHexComp,'r');
xlabel('Input Power, [MW]')
ylabel('% Loss')
xlim([1,7])
ylim([-2,15])
legend('Hex_{New} vs Hex_{OS}, Shared BoP','Hex_{Deg} vs Hex_{OS}, Shared BoP','location','Northwest')%,'Orientation','Horizontal')

subplot(2,1,2)
plot(row_DC_S2_New(:,1),loss_DCEl_S2_NewHexComp,'--b');
hold on
plot(row_DC_S2_New(:,1),loss_DCEl_S2_DegHexComp,'--r');
% title('Loss% = (qH_2_{OS Hex} - qH_2_{New/Deg Hex})/qH_2_{OS Hex} * 100')
xlabel('Input Power, [MW]')
ylabel('% Loss')
xlim([1,7])
ylim([-2,20])
legend('Hex_{New} vs Hex_{OS}, Separate BoP','Hex_{Deg} vs Hex_{OS}, Separate BoP','location','Northwest')%,'Orientation','Horizontal')


%plots for section 3.3
loss_NewHex_qfix = ((row_DC_S2_New_qlyeFix(:,end) - plot_C_S2_NewHex_qlyefix(:,end))./row_DC_S2_New_qlyeFix(:,end)).*100;
loss_OSHex_qfix = ((row_DC_S2_OS_qlyeFix(:,end) - plot_C_S2_OSHex_qlyefix(:,end))./row_DC_S2_OS_qlyeFix(:,end)).*100;
loss_DegHex_qfix = ((row_DC_S2_Deg_qlyeFix(:,end) - plot_C_S2_DegHex_qlyefix(:,end))./row_DC_S2_Deg_qlyeFix(:,end)).*100;


set(0, 'DefaultAxesFontSize', 12,  'DefaultLineLineWidth', 1.5)
fig3 = figure('NumberTitle', 'off');
hold on
subplot(3,1,1)
plot(row_DC_S2_New(:,1),loss_NewHex,'b');
hold on
plot(row_DC_S2_New(:,1),loss_NewHex_qfix,'--b');
legend('S_2F_1HX_{New}q_{lye,var} vs S_2F_2HX_{New}q_{lye,var}','S_2F_1HX_{New}q_{lye,fix} vs S_2F_2HX_{New}q_{lye,fix}','location','north')%,'Orientation','Horizontal')
xlabel('Input Power, [MW]')
ylabel('% Loss')
xlim([1,7])
ylim([-2,30])
grid on
subplot(3,1,2)
plot(row_DC_S2_New(:,1),loss_OSHex,'k');
hold on
plot(row_DC_S2_New(:,1),loss_OSHex_qfix,'--k');
legend('S_2F_1HX_{OS}q_{lye,var} vs S_2F_2HX_{OS}q_{lye,var}','S_2F_1HX_{OS}q_{lye,fix} vs S_2F_2HX_{OS}q_{lye,fix}','location','north')%,'Orientation','Horizontal')
xlabel('Input Power, [MW]')
ylabel('% Loss')
xlim([1,7])
ylim([-2,30])
grid on
subplot(3,1,3)
plot(row_DC_S2_New(:,1),loss_DegHex,'r');
hold on
plot(row_DC_S2_New(:,1),loss_DegHex_qfix,'--r');
legend('S_2F_1HX_{Deg}q_{lye,var} vs S_2F_2HX_{Deg}q_{lye,var}','S_2F_1HX_{Deg}q_{lye,fix} vs S_2F_2HX_{Deg}q_{lye,fix}','location','north')%,'Orientation','Horizontal')
xlabel('Input Power, [MW]')
ylabel('% Loss')
grid on
xlim([1,7])
ylim([-2,30])
saveas(fig3,'lyefixLoss','epsc');

%loss of changing lye flowrate from variable to fix
loss_CEl_S2_NewHex_qcomp = ((plot_C_S2_NewHex(:,end) - plot_C_S2_NewHex_qlyefix(:,end))./plot_C_S2_NewHex(:,end)).*100;
loss_DCEl_S2_NewHex_qcomp = ((row_DC_S2_New(:,end) - row_DC_S2_New_qlyeFix(:,end))./row_DC_S2_New(:,end)).*100;
loss_CEl_S2_OSHex_qcomp = ((plot_C_S2_OSHex(:,end) - plot_C_S2_OSHex_qlyefix(:,end))./plot_C_S2_OSHex(:,end)).*100;
loss_DCEl_S2_OSHex_qcomp = ((row_DC_S2_OS(:,end) - row_DC_S2_OS_qlyeFix(:,end))./row_DC_S2_OS(:,end)).*100;
loss_CEl_S2_DegHex_qcomp = ((plot_C_S2_DegHex(:,end) - plot_C_S2_DegHex_qlyefix(:,end))./plot_C_S2_DegHex(:,end)).*100;
loss_DCEl_S2_DegHex_qcomp = ((row_DC_S2_Deg(:,end) - row_DC_S2_Deg_qlyeFix(:,end))./row_DC_S2_Deg(:,end)).*100;

set(0, 'DefaultAxesFontSize', 16,  'DefaultLineLineWidth', 1.5)
fig4 = figure('NumberTitle', 'off', 'Name', 'Effect of fixing lye flowrate on the flowsheets for all the Hex designs');
hold on
subplot(3,1,1)
plot(row_DC_S2_New(:,1),loss_CEl_S2_NewHex_qcomp,'b');
hold on
plot(row_DC_S2_New(:,1),loss_DCEl_S2_NewHex_qcomp,'--b');
legend('q_{lye, fix} vs q_{lye, var}, Hex_{New}, Shared BoP','q_{lye, fix} vs q_{lye, var}, Hex_{New}, Separate BoP','location','Northwest')%,'Orientation','Horizontal')
% title('Loss = ((Varying q_{lye} - Fix q_{lye})/Varying q_{lye})*100')
xlabel('Input Power,[MW]')
ylabel('Loss %')
ylim([-2 25])
xlim([1,7])
subplot(3,1,2)
plot(row_DC_S2_New(:,1),loss_CEl_S2_OSHex_qcomp,'k');
hold on
plot(row_DC_S2_New(:,1),loss_DCEl_S2_OSHex_qcomp,'--k');
legend('q_{lye, fix} vs q_{lye, var}, Hex_{OS}, Shared BoP','q_{lye, fix} vs q_{lye, var}, Hex_{OS}, Separate BoP','location','Northwest')%,'Orientation','Horizontal')
% title('Loss = ((Varying q_{lye} - Fix q_{lye})/Varying q_{lye})*100')
xlabel('Input Power,[MW]')
ylabel('Loss %')
ylim([-2 25])
xlim([1,7])
subplot(3,1,3)
plot(row_DC_S2_New(:,1),loss_CEl_S2_DegHex_qcomp,'r');
hold on
plot(row_DC_S2_New(:,1),loss_DCEl_S2_DegHex_qcomp,'--r');
legend('q_{lye, fix} vs q_{lye, var}, Hex_{Deg}, Shared BoP','q_{lye, fix} vs q_{lye, var}, Hex_{Deg}, Separate BoP','location','Northwest')%,'Orientation','Horizontal')
% title('Loss = ((Varying q_{lye} - Fix q_{lye})/Varying q_{lye})*100')
xlabel('Input Power,[MW]')
ylabel('Loss %')
ylim([-2 25])
% saveas(fig4,'S2_C&DCEl_NewHex_qlyeComp','epsc')