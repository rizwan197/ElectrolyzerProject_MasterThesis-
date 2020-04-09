clc
clear 
close all

load Data_DecoupledElectrolyzers_State1
load Data_DecoupledElectrolyzers_State2

par.N = 3;

%active region graph for Flowsheet 2 State 1
set(0, 'DefaultAxesFontSize', 12,  'DefaultLineLineWidth', 1.5)
fig1 = figure('NumberTitle', 'off', 'Name', 'Active region map, Flowsheet 2, State 1');
hold on
linS = {'*','o','d'};
% solid line, max constraint, -- line min constraint
for i=1:par.N
    txt = ['El ',num2str(i)];
    plot(row_DC_S1(:,1),ac_DC_S1(:,i),strcat('-',linS{i},'b'),'DisplayName',strcat(txt,' Iden_M_a_x'))
    hold on
    plot(row_DC_S1(:,1),ac_DC_S1(:,par.N+i),strcat('--',linS{i},'b'),'DisplayName',strcat(txt,' Iden_M_i_n'))
    hold on
    plot(row_DC_S1(:,1),ac_DC_S1(:,2*par.N+i),strcat('-',linS{i},'r'),'DisplayName',strcat(txt,' \Delta T_M_a_x'))
    hold on
    plot(row_DC_S1(:,1),ac_DC_S1(:,6*par.N+i),strcat('--',linS{i},'k'),'DisplayName',strcat(txt,' \Delta T2_H_e_x_ _M_i_n'))
    hold on
    plot(row_DC_S1(:,1),ac_DC_S1(:,8*par.N+i),strcat('--',linS{i},'m'),'DisplayName',strcat(txt,' qlye_M_i_n'))  
end
legend show
xlabel('Supplied power, P_n_e_t [MW]')
ylabel('c/c_m_a_x')
saveas(fig1,'AC_F2S1','epsc');

%active region graph for Flowsheet 2 State 2
set(0, 'DefaultAxesFontSize', 12,  'DefaultLineLineWidth', 1.5)
fig2 = figure('NumberTitle', 'off', 'Name', 'Active region map, Flowsheet 2, State 2');
hold on
linS = {'*','o','d'};
for i=1:par.N
    txt = ['El ',num2str(i)];
    plot(row_DC_S2(:,1),ac_DC_S2(:,i),strcat('-',linS{i},'b'),'DisplayName',strcat(txt,' Iden_M_a_x'))
    hold on
    plot(row_DC_S2(:,1),ac_DC_S2(:,par.N+i),strcat('--',linS{i},'b'),'DisplayName',strcat(txt,' Iden_M_i_n'))
    hold on
    plot(row_DC_S2(:,1),ac_DC_S2(:,2*par.N+i),strcat('-',linS{i},'r'),'DisplayName',strcat(txt,' \Delta T_M_a_x'))
    hold on
    plot(row_DC_S2(:,1),ac_DC_S2(:,6*par.N+i),strcat('--',linS{i},'k'),'DisplayName',strcat(txt,' \Delta T2_H_e_x_ _M_i_n'))
    hold on
    plot(row_DC_S2(:,1),ac_DC_S2(:,7*par.N+i),strcat('-',linS{i},'m'),'DisplayName',strcat(txt,' qlye_M_a_x'))
    hold on
    plot(row_DC_S2(:,1),ac_DC_S2(:,8*par.N+i),strcat('--',linS{i},'m'),'DisplayName',strcat(txt,' qlye_M_i_n'))
    hold on
    plot(row_DC_S2(:,1),ac_DC_S2(:,9*par.N+i),strcat('-',linS{i},'g'),'DisplayName',strcat(txt,' qcw_M_a_x'))
end
legend show
xlabel('Supplied power, P_n_e_t [MW]')
ylabel('c/c_m_a_x')
saveas(fig2,'AC_F2S2','epsc');