function [z_guess,x_guess,u_guess] = init0(N,Pnet)
%function file for the initial guess

par = parElectrolyzer(N);

%algebriac state variables('z')
u_k0 = 1.8*ones(1,par.N);               %initial guess for cell voltage
P_k0 = Pnet/par.N*ones(1,par.N);        %intial guess is power divided equally among electrolyzers
i_k0 = P_k0./(u_k0.*par.EL(1).nc);      %initial guess for current
Feff_k0 = 0.97*ones(1,par.N);
nH2_k0 = 6*ones(1,par.N);               %[mol/s]
qH2Oloss_k0 = nH2_k0*par.Const.Mwt.*ones(1,par.N);%[g/s]
nH2El_net0 = sum(nH2_k0);               %[mol/s]
nH2out_net0 = sum(nH2_k0);
nO2El_net0 = 0.5*nH2El_net0;
nO2out_net0 = 0.5*nH2out_net0;
T_El_out0 = 80;                         %initial guess for the temperature of lye entering in the heat exchanger

%differential state variables('x')
T_k0 = 75*ones(1,par.N);
Psto_H20 = 25;        %initial H2 storage pressure (calculated from steady state solution) [bar]
Psto_O20 = 25;        %initial O2 storage pressure (calculated from steady state solution) [bar]
Mass_Bt0 = 6000000;  %mass of the liquid in the buffer tank,[g] 6000kg
T_bt_out0 = 70;       %Initial guess for the temperature of lye mixture at the exit of the buffer tank,[degC]
T_El_in0 = 65;        %initial guess for the temperature of inlet lye into the electrolyzer, [deg C]
T_cw_out0 = 20;       %initial guess for the exit temperature of the cooling water leaving heat exchanger,[deg C]

z_guess = [u_k0 i_k0 P_k0 Feff_k0 nH2_k0 qH2Oloss_k0 nH2El_net0 nH2out_net0 nO2El_net0 nO2out_net0 T_El_out0];
x_guess = [T_k0 Psto_H20 Psto_O20 Mass_Bt0 T_bt_out0 T_El_in0 T_cw_out0];

%initial guess for input variables('u')
U_El_k_0 = 414.0301;      %voltage across electrolyzers, [Volts]
q_lye_k_0 = 6648*ones(1,par.N);         %lye flowrate, [g/s]
q_cw_0 = 2.0698e4;                      %cooling water flow rate, [g/s]
zH2_0 = 0.4;
zO2_0 = 0.4;
q_H2O_0 = 324.2657;                     %total water lost during electrolysis, [grams/sec]

u_guess = [U_El_k_0 q_lye_k_0 q_cw_0 zH2_0 zO2_0 q_H2O_0];

end
