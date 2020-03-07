clc
clear
close all

load data_qlye2step3000_MVQcool_12hr.mat

%% Plotting the results
figure()
plot(I)
xlabel('Time, s')
ylabel('Current, A')
legend('El 1','El 2', 'El 3')
grid on

figure()
plot(U)
xlabel('Time, s')
ylabel('Cell voltage, V/cell')
legend('El 1','El 2', 'El 3')
grid on

figure()
subplot(1,4,1)
plot(P(:,1),'b')
xlabel('Time, s')
ylabel('P_1, Watts')
grid on
subplot(1,4,2)
plot(P(:,2),'k')
xlabel('Time, s')
ylabel('P_2, Watts')
grid on
subplot(1,4,3)
plot(P(:,3),'r')
xlabel('Time, s')
ylabel('P_3, Watts')
grid on
subplot(1,4,4)
plot(P_net)
xlabel('Time, s')
ylabel('Power, Watts')
grid on

figure()
plot(nH2elout)
hold on
plot(nH2in)
xlabel('Time, s')
ylabel('H_2 production rate, mol/s')
legend('El 1','El 2', 'El 3','H_2_n_e_t El')
grid on

figure()
subplot(3,1,1)
plot(Temp(2:end,1))
xlabel('Time, s')
ylabel('T_1, C')
grid on

subplot(3,1,2)
plot(Temp(2:end,2))
xlabel('Time, s')
ylabel('T_2, C')
grid on

subplot(3,1,3)
plot(Temp(2:end,3))
xlabel('Time, s')
ylabel('T_3, C')
grid on

figure()

subplot(2,1,1)
plot(Tout)
xlabel('Time, s')
ylabel('T_o_u_t, C')
grid on

subplot(2,1,2)
plot(Tin)
xlabel('Time, s')
ylabel('T_i_n, C')
grid on

figure()
subplot(1,2,1)
plot(PstoH2)
xlabel('t, s')
ylabel('Psto_H_2 bar')
grid on

subplot(1,2,2)
plot(PstoO2)
xlabel('t, s')
ylabel('Psto_O_2 bar')
grid on

%plot of the step in MV

figure()
%plot(V_El)
plot(Qc)
xlabel('Time, s')
ylabel('Q_c_o_o_l, Watts')
grid on
%legend('El 1','El 2', 'El 3')