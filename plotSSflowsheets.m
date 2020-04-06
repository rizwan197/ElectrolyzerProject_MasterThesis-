clc
clear 
close all

load Data_CoupledElectrolyzer

set(0, 'DefaultAxesFontSize', 12,  'DefaultLineLineWidth', 1.5)
fig1 = figure('NumberTitle', 'off', 'Name', 'Electrolyzer current');
hold on
plot(row1(:,2),row1(:,end))
grid on
saveas(fig1,'FlowsheetResult','epsc');
