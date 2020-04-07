clc
clear 
close all

load Data_CoupledElectrolyzerState2
load Data_DecoupledElectrolyzers_State2
par.N = 3;
% set(0, 'DefaultAxesFontSize', 12,  'DefaultLineLineWidth', 1.5)
% fig1 = figure('NumberTitle', 'off', 'Name', 'Electrolyzer current');
% hold on
% plot(row_C_S2(:,2),row_C_S2(:,end))%plot for lower Pnet needs to be fixed
% hold on
% plot(row_DC_S2(:,2),row_DC_S2(:,end))
% grid on
% legend('Flowsheet 1','Flowsheet 2')
% saveas(fig1,'FlowsheetResult','epsc');

%active region graph for Flowsheet 1 State 2
set(0, 'DefaultAxesFontSize', 12,  'DefaultLineLineWidth', 1.5)
fig1 = figure('NumberTitle', 'off', 'Name', 'Electrolyzer current');
hold on

plot(row_C_S2(:,1),ac_C_S2(:,[1,par.N+1,2*par.N+1,3*par.N+1,5*par.N+2,5*par.N+3,6*par.N+3,7*par.N+3]))
legend('Iden_M_a_x','Iden_M_i_n','T1_M_a_x','T1_M_i_n','\delta T2Hex_M_i_n','qlye1_M_a_x','qlye1_M_i_n','qcw_M_a_x')

saveas(fig1,'AC_F1S2','epsc');
