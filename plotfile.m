close all
load data_Step_qcw_40hr

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
subplot(2,1,1)
plot(V_H2)
xlabel('Time, s')
ylabel('H_2 production rate, Nm^3/hr')
legend('El 1','El 2', 'El 3')
grid on

subplot(2,1,2)
plot(Ps)
xlabel('Time, s')
ylabel('Specific electricity consumption, kWh/Nm^3')
legend('El 1','El 2', 'El 3')
grid on

figure()
subplot(2,1,1)
plot(PstoH2)
xlabel('Time, s')
ylabel('H_2 Storage pressure, bar')
grid on
subplot(2,1,2)
plot(PcompH2)
xlabel('Time, s')
ylabel('H_2 compressor power, watts')
grid on

figure()
subplot(2,1,1)
plot(Temp(:,1))
hold on
plot(Temp(:,2))
hold on
plot(Temp(:,3))
xlabel('Time, s')
ylabel('T_k, C')
grid on
subplot(2,1,2)
plot(Telout)
xlabel('Time, s')
ylabel('T_o_u_t, C')
grid on

figure()
subplot(3,1,1)
plot(qlye(:,1))
xlabel('Time, s')
ylabel('Lye flowrate to El 1')
grid on
subplot(3,1,2)
plot(q_cw)
xlabel('Time, s')
ylabel('Cooling duty')
grid on
subplot(3,1,3)
plot(Telin,'k')
hold on
plot(T_El_in_set,'r--')
xlabel('Time, s')
ylabel('T_E_l_ _i_n, C')
%ylim([63.5 65.5])
grid on