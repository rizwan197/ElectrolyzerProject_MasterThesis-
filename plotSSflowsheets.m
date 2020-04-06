clc
clear 
close all

load Data_CoupledElectrolyzer
load Data_DecoupledElectrolyzers

set(0, 'DefaultAxesFontSize', 12,  'DefaultLineLineWidth', 1.5)
fig1 = figure('NumberTitle', 'off', 'Name', 'Electrolyzer current');
hold on
plot(row1(:,2),row1(:,end))%plot for lower Pnet needs to be fixed
hold on
plot(row(:,2),row(:,end))
grid on
legend('Flowsheet 1','Flowsheet 2')
saveas(fig1,'FlowsheetResult','epsc');
