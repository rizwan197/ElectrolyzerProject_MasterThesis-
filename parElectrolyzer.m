function par = parElectrolyzer(N)

%This script defines values of the input parameters for all electrolyzers.

par.Const = struct('ze',2,'FC',96485,'R',8.314,'Cp',4.186,'CpLye',3.1006,...
    'Mwt',18,'MwtH2',2.01588,'Tref',25,'rho',1000,'rhoLye',1258.2,'Vc',2.0681,'Vh',1.9944);
%Cp=specific heat of water, [J/gK];Mwt=mol. wt of H2O, rho=density of
%water/lye[kg/m3],Vc=volume of cold side of heat exchanger[m3],Vh=volume of
%hot side of heat exchanger[m3]

par.Comp = struct('alpha',0.63,'k',1.62,'Tel',25+273,'Pel',3);
par.Storage = struct('VstoH2',965000,'VstoO2',482500,'PoutH2',19,'PoutO2',19,...
    'TstoH2',25+273.15,'TstoO2',25+273.15,'Rg',8.314e-2,'VdispH2',0.5,'VdispO2',0.5);
%VstoH2 and VstoO2 are in litres

par.Tw_in = 10;             %inlet temperature of the cooling water in lye circulation heat exchanger
par.Hex.UA = 20.48e3;%1*1.5205e4; %UA of heat exchanger [W/K], calculated from previous ss formulation
par.kvalveH2 = 14.723;      %valve constant for the outlet valve of hydrogen storage tank, calculated for 25 bar storage pressure at SS   
par.kvalveO2 = 7.362;       %valve constant for the outlet valve of oxygen storage tank, calculated for 25 bar storage pressure at SS    
par.sigma = 5.672*10^-8;    %stefan-boltzmann constant [W/m^2 K^4]
par.em = 0.8;               %emissivity [-]

%% Parameters for U-I relationship in Ulleberg's model
par.U = struct([]);
par.TherMo = struct([]);
par.EL = struct([]);
for i =1:N
    %Vidar's Parameters
    par.U(i).r1 = 0.000218155;         %ohm m^2
    par.U(i).r2 = -0.000000425;          %ohm m^2 C^-1
    par.U(i).s = 0.1179375;              %Vs
    par.U(i).t1 = -0.14529;              %A^-1 m^2
    par.U(i).t2 = 11.794;                %A^-1 m^2 C^-1
    par.U(i).t3 = 395.68;                %A^-1 m^2 C^-2
    par.U(i).f1 = 120;                   %mA^2 cm^-4
    par.U(i).f2 = 0.98;                  %dimensionless
    
    %Ulleberg parameters
%     par.U(i).r1 = 0.0000805;
%     par.U(i).r2 = -0.00000025;
%     par.U(i).s = 0.185;
%     par.U(i).t1 = -0.1002;
%     par.U(i).t2 = 8.424;
%     par.U(i).t3 = 247.3;
%     par.U(i).f1 = 250;
%     par.U(i).f2 = 0.98;
    
    %% Parameters for the thermal model
%     par.TherMo(i).CtS = 625/27;         %Suggested by Vidar, specific thermal capacity of electrolyzer i.e. Ct/P, [kJ/kWatts*C] 
    par.TherMo(i).Ct = 625/27*2134;     %Cts*Pnom [kJ/s], Assuming Pnom = 2134 kWatts                   
    par.TherMo(i).hc = 5.5;             %convective heat transfer coefficient W/m^2 C
%     par.TherMo(i).A_surf = 0.1;         %specific radition area per kA current per cell, [m^2/kA*Ncell]
    par.TherMo(i).A_El = 0.1*5.72*230;  %surface area of the electrolyzer, A_surf*Inom*Ncell [m^2]
    
    %% Parameters for Faraday effeciency calculations
    par.EL(i).Utn = 1.482;             %thermoneutral voltage, [V]
    par.EL(i).nc = 230;                %no. of cells
    par.EL(i).A = 2.6;                 %electrode area of each cell, [m^2]
    par.EL(i).Ta = 20;                 %ambient temp, [C]
    par.EL(i).Tstd = 25;               %standard temperature, [C]
    
end
% El #2, performing at 85% of electrolyzer 1 
par.U(2).r1 = par.U(2).r1*1.2;             %ohm m^2
par.U(2).s = par.U(2).s*1.2;               %V
par.U(2).f1 = par.U(2).f1*1.2;             %mA^2 cm^-4
par.U(2).f2 = 0.97;

%El #3, performing at 70% of electrolyzer 1
par.U(3).r1 = par.U(3).r1*1.3;              %ohm m^2
par.U(3).s = par.U(3).s*1.3;                %V
par.U(3).f1 = par.U(3).f1*1.3;              %mA^2 cm^-4
par.U(3).f2 = 0.96;

par.N=N;
end