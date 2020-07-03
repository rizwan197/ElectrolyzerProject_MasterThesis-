clc
clear
close all

% load Data_CoupledElectrolyzerState2
% load Data_CoupledElectrolyzerState1
% load Data_DecoupledElectrolyzers_State1
% load Data_DecoupledElectrolyzers_State2

load DataThesis_ActiveRegionPlot

par.N = 3;

%% active region graphs for the thesis
%  ac_C_S2(counter,:) = [Iden/198.5, 32./Iden, Tk/80, 25./Tk, (max(Tk)-T_El_in_set)/30, Pcons/Pnet,...
%         qlye_kgs/10, 0.5./qlye_kgs, qcw_kgs/80 1e-5/qcw_kgs];

set(0, 'DefaultAxesFontSize', 16,  'DefaultLineLineWidth', 4)
fig1 = figure('NumberTitle', 'off', 'Name', 'Active region map, Flowsheet 1, State 2');
hold on
linS = {'*','o','d'};
for i=1:par.N
    txt = [',',num2str(i)];
    plot(row_C_S2(:,1),ac_C_S2(:,i),strcat('-',linS{i},'b'),'DisplayName',strcat('g_{',num2str(5*par.N+i),'} (',' I_{den',txt,'}/I_{den',txt,'}^{max})'), 'MarkerSize',15);%Iden_max constraint
    hold on
end
for i=1:par.N
    txt = num2str(i);
    plot(row_C_S2(:,1),ac_C_S2(:,2*par.N+i),strcat('-',linS{i},'r'),'DisplayName',strcat('g_{',num2str(7*par.N+i),'} (', ' T_{',txt,'}/T_{',txt,'}^{max})'), 'MarkerSize',15);%T_max constraint
    hold on
end
hold on
plot(row_C_S2(:,1),ac_C_S2(:,4*par.N+1),strcat('-','k'),'DisplayName',strcat('g_{25} (','\Delta T_{El,in}/ \Delta T_{El,in}^{max})'))%Tk-Telin constraint
hold on
plot(row_C_S2(:,1),ac_C_S2(:,4*par.N+2),strcat('-','g'),'DisplayName',strcat('g_{26} (','P_{cons}/P_{net})'))%Pcons constraint

ylim([0 1.2])
xlim([3.4,7])
legend show
legend('FontSize',18)
% legend('Orientation','horizontal')
xlabel('Supplied power, P_n_e_t [MW]')
ylabel('Scaled output (CV) constraint, g/g_{max}')
h1 = xline(5.1,'--k','Medium Power','LabelOrientation','Horizontal','FontSize',18);
h1.Annotation.LegendInformation.IconDisplayStyle = 'off';
h2 = xline(5.8,'--k','High Power','LabelOrientation','Horizontal','FontSize',18);
h2.Annotation.LegendInformation.IconDisplayStyle = 'off';
h3 = xline(3.4,'--k','Low Power','LabelOrientation','Horizontal','FontSize',18);
h3.Annotation.LegendInformation.IconDisplayStyle = 'off';
box on
saveas(fig1,'ThesisAC_Output','epsc');

set(0, 'DefaultAxesFontSize', 16,  'DefaultLineLineWidth', 4)
fig2 = figure('NumberTitle', 'off', 'Name', 'Active region map, Flowsheet 1, State 2');
hold on
linS = {'*','o','d'};
for i=1:par.N
    txt = [',',num2str(i)];
    plot(row_C_S2(:,1),ac_C_S2(:,4*par.N+2+i),strcat('-',linS{i},'m'),'DisplayName',strcat('g_{',num2str(par.N+i),'} (',' q_{lye',txt,'}/q_{lye',txt,'}^{max})'), 'MarkerSize',15);%qlye_max constraint
    hold on
end
ylim([0 1.2])
xlim([3.4,7])
legend show
legend('FontSize',18)
% legend('Orientation','horizontal')
xlabel('Supplied power, P_n_e_t [MW]')
ylabel('Scaled input (MV) constraint, g/g_{max}')
h4 = xline(5.1,'--k','Medium Power','DisplayName','P_{net} = 5.1 MW','LabelOrientation','Horizontal','FontSize',18);
h4.Annotation.LegendInformation.IconDisplayStyle = 'off';
h5 = xline(5.8,'--k','High Power','DisplayName','P_{net} = 5.8 MW','LabelOrientation','Horizontal','FontSize',18);
h5.Annotation.LegendInformation.IconDisplayStyle = 'off';
h6 = xline(3.4,'--k','Low Power','DisplayName','P_{net} = 3.4 MW','LabelOrientation','Horizontal','FontSize',18);
h6.Annotation.LegendInformation.IconDisplayStyle = 'off';
box on
saveas(fig2,'ThesisAC_Input','epsc');

%% %active region graph for Flowsheet 1 State 1
% set(0, 'DefaultAxesFontSize', 12,  'DefaultLineLineWidth', 1.5)
% fig1 = figure('NumberTitle', 'off', 'Name', 'Active region map, Flowsheet 1, State 1');
% hold on
% subplot(2,1,1)
% linS = {'*','o','d'};
% for i=1:par.N
%     txt = ['El ',num2str(i)];
%     plot(row_C_S1(:,1),ac_C_S1(:,i),strcat('-',linS{i},'b'),'DisplayName',strcat(txt,' Iden_M_a_x'))
%     hold on
%     plot(row_C_S1(:,1),ac_C_S1(:,par.N+i),strcat('--',linS{i},'b'),'DisplayName',strcat(txt,' Iden_M_i_n'))
%     hold on
%     plot(row_C_S1(:,1),ac_C_S1(:,2*par.N+i),strcat('-',linS{i},'r'),'DisplayName',strcat(txt,' T_M_a_x'))
% end
% hold on
% plot(row_C_S1(:,1),ac_C_S1(:,5*par.N+2),strcat('--','k'),'DisplayName',' \Delta T2_H_e_x_ _M_i_n')
% ylim([0 1.2])
% xlim([0,7])
% legend show
% legend('Orientation','horizontal')
% xlabel('Supplied power, P_n_e_t [MW]')
% ylabel('Scaled active constraint, c/c_m_a_x')
% xline(1,'-k','LineWidth',1.2,'DisplayName','P_n_e_t_ _M_i_n');
% 
% subplot(2,1,2)
% linS = {'*','o','d'};
% for i=1:par.N
%     txt = ['El ',num2str(i)];
%     plot(row_C_S1(:,1),ac_C_S1(:,6*par.N+2+i),strcat('--',linS{i},'m'),'DisplayName',strcat(txt,' qlye_M_i_n'))
%     hold on
% end
% ylim([0 1.2])
% xlim([0,7])
% legend show
% legend('Orientation','horizontal')
% xlabel('Supplied power, P_n_e_t [MW]')
% ylabel('Scaled active constraint, c/c_m_a_x')
% xline(1,'-k','LineWidth',1.2,'DisplayName','P_n_e_t_ _M_i_n');
% saveas(fig1,'AC_F1S1','epsc');

%% %active region graph for Flowsheet 1 State 2
% set(0, 'DefaultAxesFontSize', 12,  'DefaultLineLineWidth', 1.5)
% fig2 = figure('NumberTitle', 'off', 'Name', 'Active region map, Flowsheet 1, State 2');
% hold on
% subplot(2,1,1)
% linS = {'*','o','d'};
% for i=1:par.N
%     txt = ['El ',num2str(i)];
%     plot(row_C_S2(:,1),ac_C_S2(:,i),strcat('-',linS{i},'b'),'DisplayName',strcat(txt,' Iden_M_a_x'))
%     hold on
%     plot(row_C_S2(:,1),ac_C_S2(:,par.N+i),strcat('--',linS{i},'b'),'DisplayName',strcat(txt,' Iden_M_i_n'))
%     hold on
%     plot(row_C_S2(:,1),ac_C_S2(:,2*par.N+i),strcat('-',linS{i},'r'),'DisplayName',strcat(txt,' T_M_a_x'))
%     hold on
% end
% hold on
% plot(row_C_S2(:,1),ac_C_S2(:,5*par.N+2),strcat('--','k'),'DisplayName',' \Delta T2_H_e_x_ _M_i_n')
% ylim([0 1.2])
% xlim([0,7])
% legend show
% legend('Orientation','horizontal')
% xlabel('Supplied power, P_n_e_t [MW]')
% ylabel('Scaled active constraint, c/c_m_a_x')
% xline(1.2,'-k','LineWidth',1.2,'DisplayName','P_n_e_t_ _M_i_n');
% subplot(2,1,2)
% linS = {'*','o','d'};
% for i=1:par.N
%     txt = ['El ',num2str(i)];
%     plot(row_C_S2(:,1),ac_C_S2(:,5*par.N+2+i),strcat('-',linS{i},'m'),'DisplayName',strcat(txt,' qlye_M_a_x'))
%     hold on
%     plot(row_C_S2(:,1),ac_C_S2(:,6*par.N+2+i),strcat('--',linS{i},'m'),'DisplayName',strcat(txt,' qlye_M_i_n'))
%     hold on
% end
% plot(row_C_S2(:,1),ac_C_S2(:,7*par.N+3),strcat('-','g'),'DisplayName',' qcw_M_a_x')
% ylim([0 1.2])
% xlim([0,7])
% legend show
% legend('Orientation','horizontal')
% xlabel('Supplied power, P_n_e_t [MW]')
% ylabel('Scaled active constraint, c/c_m_a_x')
% xline(1.2,'-k','LineWidth',1.2,'DisplayName','P_n_e_t_ _M_i_n');
% saveas(fig2,'AC_F1S2','epsc');

% %% active region graph for Flowsheet 2 State 1
% set(0, 'DefaultAxesFontSize', 12,  'DefaultLineLineWidth', 1.5)
% fig3 = figure('NumberTitle', 'off', 'Name', 'Active region map, Flowsheet 2, State 1');
% hold on
% subplot(2,1,1)
% linS = {'*','o','d'};
% % solid line, max constraint, -- line min constraint
% for i=1:par.N
%     txt = ['El ',num2str(i)];
%     plot(row_DC_S1(:,1),ac_DC_S1(:,i),strcat('-',linS{i},'b'),'DisplayName',strcat(txt,' Iden_M_a_x'))
%     hold on
%     plot(row_DC_S1(:,1),ac_DC_S1(:,par.N+i),strcat('--',linS{i},'b'),'DisplayName',strcat(txt,' Iden_M_i_n'))
%     hold on
%     plot(row_DC_S1(:,1),ac_DC_S1(:,2*par.N+i),strcat('-',linS{i},'r'),'DisplayName',strcat(txt,' T_M_a_x'))
%     hold on
%     plot(row_DC_S1(:,1),ac_DC_S1(:,6*par.N+i),strcat('--',linS{i},'k'),'DisplayName',strcat(txt,' \Delta T2_H_e_x_ _M_i_n'))
%     hold on
% end
% ylim([0 1.2])
% xlim([0,7])
% legend show
% legend('Orientation','horizontal')
% xlabel('Supplied power, P_n_e_t [MW]')
% ylabel('Scaled active constraint, c/c_m_a_x')
% xline(1,'-k','LineWidth',1.2,'DisplayName','P_n_e_t_ _M_i_n');
% 
% subplot(2,1,2)
% linS = {'*','o','d'};
% % solid line, max constraint, -- line min constraint
% for i=1:par.N
%     txt = ['El ',num2str(i)];
%     plot(row_DC_S1(:,1),ac_DC_S1(:,8*par.N+i),strcat('--',linS{i},'m'),'DisplayName',strcat(txt,' qlye_M_i_n'))  
%     hold on
% end
% ylim([0 1.2])
% xlim([0,7])
% legend show
% legend('Orientation','horizontal')
% xlabel('Supplied power, P_n_e_t [MW]')
% ylabel('Scaled active constraint, c/c_m_a_x')
% xline(1,'-k','LineWidth',1.2,'DisplayName','P_n_e_t_ _M_i_n');
% saveas(fig3,'AC_F2S1','epsc');
% 
% %active region graph for Flowsheet 2 State 2
% set(0, 'DefaultAxesFontSize', 12,  'DefaultLineLineWidth', 1.5)
% fig4 = figure('NumberTitle', 'off', 'Name', 'Active region map, Flowsheet 2, State 2');
% hold on
% subplot(2,1,1)
% linS = {'*','o','d'};
% for i=1:par.N
%     txt = ['El ',num2str(i)];
%     plot(row_DC_S2(:,1),ac_DC_S2(:,i),strcat('-',linS{i},'b'),'DisplayName',strcat(txt,' Iden_M_a_x'))
%     hold on
%     plot(row_DC_S2(:,1),ac_DC_S2(:,par.N+i),strcat('--',linS{i},'b'),'DisplayName',strcat(txt,' Iden_M_i_n'))
%     hold on
%     plot(row_DC_S2(:,1),ac_DC_S2(:,2*par.N+i),strcat('-',linS{i},'r'),'DisplayName',strcat(txt,' T_M_a_x'))
%     hold on
%     plot(row_DC_S2(:,1),ac_DC_S2(:,6*par.N+i),strcat('--',linS{i},'k'),'DisplayName',strcat(txt,' \Delta T2_H_e_x_ _M_i_n'))
%     hold on
% end
% ylim([0 1.2])
% xlim([0,7])
% legend show
% legend('Orientation','horizontal')
% xlabel('Supplied power, P_n_e_t [MW]')
% ylabel('Scaled active constraint, c/c_m_a_x')
% xline(1.1,'-k','LineWidth',1.2,'DisplayName','P_n_e_t_ _M_i_n');
% 
% subplot(2,1,2)
% linS = {'*','o','d'};
% for i=1:par.N
%     txt = ['El ',num2str(i)];
%     plot(row_DC_S2(:,1),ac_DC_S2(:,7*par.N+i),strcat('-',linS{i},'m'),'DisplayName',strcat(txt,' qlye_M_a_x'))
%     hold on
%     plot(row_DC_S2(:,1),ac_DC_S2(:,8*par.N+i),strcat('--',linS{i},'m'),'DisplayName',strcat(txt,' qlye_M_i_n'))
%     hold on
%     plot(row_DC_S2(:,1),ac_DC_S2(:,9*par.N+i),strcat('-',linS{i},'g'),'DisplayName',strcat(txt,' qcw_M_a_x'))
%     hold on
% end
% ylim([0 1.2])
% xlim([0,7])
% legend show
% legend('Orientation','horizontal')
% xlabel('Supplied power, P_n_e_t [MW]')
% ylabel('Scaled active constraint, c/c_m_a_x')
% xline(1.1,'-k','LineWidth',1.2,'DisplayName','P_n_e_t_ _M_i_n');
% saveas(fig4,'AC_F2S2','epsc');