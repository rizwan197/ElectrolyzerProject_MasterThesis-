clc
clear
close all

load Data_CoupledFlowsheet
load Data_DecoupledFlowsheet

set(0, 'DefaultAxesFontSize', 12,  'DefaultLineLineWidth', 1.5)
fig1 = figure('NumberTitle', 'off', 'Name', 'Electrolyzer current');
hold on
vH2coupled = sum(row1(16:18));
plot(row1(:,2),vH2coupled)
hold on
plot(Power,Flowsheet2)
xlabel('P_n_e_t, Megawatts')
ylabel('H_2 production, Nm^3/hr')
legend('Flowsheet 1', 'Flowsheet 2','Location','northwest')
saveas(fig1,'FlowsheetResults','epsc');