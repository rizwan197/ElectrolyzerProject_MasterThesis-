clc
clear
close all

load Data_CoupledElectrolyzerState2
load Data_CoupledElectrolyzerState1
load Data_DecoupledElectrolyzers_State1
load Data_DecoupledElectrolyzers_State2

par.N = 3;

%active region graph for Flowsheet 1 State 1
set(0, 'DefaultAxesFontSize', 12,  'DefaultLineLineWidth', 1.5)
fig1 = figure('NumberTitle', 'off', 'Name', 'Active region map, Flowsheet 1, State 1');
hold on
subplot(2,1,1)
linS = {'*','o','d'};
for i=1:par.N
    txt = ['El ',num2str(i)];
    plot(row_C_S1(:,1),ac_C_S1(:,i),strcat('-',linS{i},'b'),'DisplayName',strcat(txt,' Iden_M_a_x'))
    hold on
    plot(row_C_S1(:,1),ac_C_S1(:,par.N+i),strcat('--',linS{i},'b'),'DisplayName',strcat(txt,' Iden_M_i_n'))
    hold on
    plot(row_C_S1(:,1),ac_C_S1(:,2*par.N+i),strcat('-',linS{i},'r'),'DisplayName',strcat(txt,' T_M_a_x'))
end
hold on
plot(row_C_S1(:,1),ac_C_S1(:,5*par.N+2),strcat('--','k'),'DisplayName',' \Delta T2_H_e_x_ _M_i_n')
ylim([0 1.2])
xlim([0,7])
legend show
legend('Orientation','horizontal')
xlabel('Supplied power, P_n_e_t [MW]')
ylabel('Scaled active constraint, c/c_m_a_x')
xline(1,'-k','LineWidth',1.2,'DisplayName','P_n_e_t_ _M_i_n');

subplot(2,1,2)
linS = {'*','o','d'};
for i=1:par.N
    txt = ['El ',num2str(i)];
    plot(row_C_S1(:,1),ac_C_S1(:,6*par.N+2+i),strcat('--',linS{i},'m'),'DisplayName',strcat(txt,' qlye_M_i_n'))
    hold on
end
ylim([0 1.2])
xlim([0,7])
legend show
legend('Orientation','horizontal')
xlabel('Supplied power, P_n_e_t [MW]')
ylabel('Scaled active constraint, c/c_m_a_x')
xline(1,'-k','LineWidth',1.2,'DisplayName','P_n_e_t_ _M_i_n');
saveas(fig1,'AC_F1S1','epsc');

%active region graph for Flowsheet 1 State 2
set(0, 'DefaultAxesFontSize', 12,  'DefaultLineLineWidth', 1.5)
fig2 = figure('NumberTitle', 'off', 'Name', 'Active region map, Flowsheet 1, State 2');
hold on
subplot(2,1,1)
linS = {'*','o','d'};
for i=1:par.N
    txt = ['El ',num2str(i)];
    plot(row_C_S2(:,1),ac_C_S2(:,i),strcat('-',linS{i},'b'),'DisplayName',strcat(txt,' Iden_M_a_x'))
    hold on
    plot(row_C_S2(:,1),ac_C_S2(:,par.N+i),strcat('--',linS{i},'b'),'DisplayName',strcat(txt,' Iden_M_i_n'))
    hold on
    plot(row_C_S2(:,1),ac_C_S2(:,2*par.N+i),strcat('-',linS{i},'r'),'DisplayName',strcat(txt,' T_M_a_x'))
    hold on
end
hold on
plot(row_C_S2(:,1),ac_C_S2(:,5*par.N+2),strcat('--','k'),'DisplayName',' \Delta T2_H_e_x_ _M_i_n')
ylim([0 1.2])
xlim([0,7])
legend show
legend('Orientation','horizontal')
xlabel('Supplied power, P_n_e_t [MW]')
ylabel('Scaled active constraint, c/c_m_a_x')
xline(1.2,'-k','LineWidth',1.2,'DisplayName','P_n_e_t_ _M_i_n');
subplot(2,1,2)
linS = {'*','o','d'};
for i=1:par.N
    txt = ['El ',num2str(i)];
    plot(row_C_S2(:,1),ac_C_S2(:,5*par.N+2+i),strcat('-',linS{i},'m'),'DisplayName',strcat(txt,' qlye_M_a_x'))
    hold on
    plot(row_C_S2(:,1),ac_C_S2(:,6*par.N+2+i),strcat('--',linS{i},'m'),'DisplayName',strcat(txt,' qlye_M_i_n'))
    hold on
end
plot(row_C_S2(:,1),ac_C_S2(:,7*par.N+3),strcat('-','g'),'DisplayName',' qcw_M_a_x')
ylim([0 1.2])
xlim([0,7])
legend show
legend('Orientation','horizontal')
xlabel('Supplied power, P_n_e_t [MW]')
ylabel('Scaled active constraint, c/c_m_a_x')
xline(1.2,'-k','LineWidth',1.2,'DisplayName','P_n_e_t_ _M_i_n');
saveas(fig2,'AC_F1S2','epsc');

%% active region graph for Flowsheet 2 State 1
set(0, 'DefaultAxesFontSize', 12,  'DefaultLineLineWidth', 1.5)
fig3 = figure('NumberTitle', 'off', 'Name', 'Active region map, Flowsheet 2, State 1');
hold on
subplot(2,1,1)
linS = {'*','o','d'};
% solid line, max constraint, -- line min constraint
for i=1:par.N
    txt = ['El ',num2str(i)];
    plot(row_DC_S1(:,1),ac_DC_S1(:,i),strcat('-',linS{i},'b'),'DisplayName',strcat(txt,' Iden_M_a_x'))
    hold on
    plot(row_DC_S1(:,1),ac_DC_S1(:,par.N+i),strcat('--',linS{i},'b'),'DisplayName',strcat(txt,' Iden_M_i_n'))
    hold on
    plot(row_DC_S1(:,1),ac_DC_S1(:,2*par.N+i),strcat('-',linS{i},'r'),'DisplayName',strcat(txt,' T_M_a_x'))
    hold on
    plot(row_DC_S1(:,1),ac_DC_S1(:,6*par.N+i),strcat('--',linS{i},'k'),'DisplayName',strcat(txt,' \Delta T2_H_e_x_ _M_i_n'))
    hold on
end
ylim([0 1.2])
xlim([0,7])
legend show
legend('Orientation','horizontal')
xlabel('Supplied power, P_n_e_t [MW]')
ylabel('Scaled active constraint, c/c_m_a_x')
xline(1,'-k','LineWidth',1.2,'DisplayName','P_n_e_t_ _M_i_n');

subplot(2,1,2)
linS = {'*','o','d'};
% solid line, max constraint, -- line min constraint
for i=1:par.N
    txt = ['El ',num2str(i)];
    plot(row_DC_S1(:,1),ac_DC_S1(:,8*par.N+i),strcat('--',linS{i},'m'),'DisplayName',strcat(txt,' qlye_M_i_n'))  
    hold on
end
ylim([0 1.2])
xlim([0,7])
legend show
legend('Orientation','horizontal')
xlabel('Supplied power, P_n_e_t [MW]')
ylabel('Scaled active constraint, c/c_m_a_x')
xline(1,'-k','LineWidth',1.2,'DisplayName','P_n_e_t_ _M_i_n');
saveas(fig3,'AC_F2S1','epsc');

%active region graph for Flowsheet 2 State 2
set(0, 'DefaultAxesFontSize', 12,  'DefaultLineLineWidth', 1.5)
fig4 = figure('NumberTitle', 'off', 'Name', 'Active region map, Flowsheet 2, State 2');
hold on
subplot(2,1,1)
linS = {'*','o','d'};
for i=1:par.N
    txt = ['El ',num2str(i)];
    plot(row_DC_S2(:,1),ac_DC_S2(:,i),strcat('-',linS{i},'b'),'DisplayName',strcat(txt,' Iden_M_a_x'))
    hold on
    plot(row_DC_S2(:,1),ac_DC_S2(:,par.N+i),strcat('--',linS{i},'b'),'DisplayName',strcat(txt,' Iden_M_i_n'))
    hold on
    plot(row_DC_S2(:,1),ac_DC_S2(:,2*par.N+i),strcat('-',linS{i},'r'),'DisplayName',strcat(txt,' T_M_a_x'))
    hold on
    plot(row_DC_S2(:,1),ac_DC_S2(:,6*par.N+i),strcat('--',linS{i},'k'),'DisplayName',strcat(txt,' \Delta T2_H_e_x_ _M_i_n'))
    hold on
end
ylim([0 1.2])
xlim([0,7])
legend show
legend('Orientation','horizontal')
xlabel('Supplied power, P_n_e_t [MW]')
ylabel('Scaled active constraint, c/c_m_a_x')
xline(1.1,'-k','LineWidth',1.2,'DisplayName','P_n_e_t_ _M_i_n');

subplot(2,1,2)
linS = {'*','o','d'};
for i=1:par.N
    txt = ['El ',num2str(i)];
    plot(row_DC_S2(:,1),ac_DC_S2(:,7*par.N+i),strcat('-',linS{i},'m'),'DisplayName',strcat(txt,' qlye_M_a_x'))
    hold on
    plot(row_DC_S2(:,1),ac_DC_S2(:,8*par.N+i),strcat('--',linS{i},'m'),'DisplayName',strcat(txt,' qlye_M_i_n'))
    hold on
    plot(row_DC_S2(:,1),ac_DC_S2(:,9*par.N+i),strcat('-',linS{i},'g'),'DisplayName',strcat(txt,' qcw_M_a_x'))
    hold on
end
ylim([0 1.2])
xlim([0,7])
legend show
legend('Orientation','horizontal')
xlabel('Supplied power, P_n_e_t [MW]')
ylabel('Scaled active constraint, c/c_m_a_x')
xline(1.1,'-k','LineWidth',1.2,'DisplayName','P_n_e_t_ _M_i_n');
saveas(fig4,'AC_F2S2','epsc');