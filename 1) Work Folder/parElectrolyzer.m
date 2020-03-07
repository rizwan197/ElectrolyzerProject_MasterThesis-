function par = parElectrolyzer(N)

%This script defines values of the input parameters for all electrolyzers.

par.Const = struct('ze',2,'FC',96485,'R',8.314,'Cp',4.186,'CpLye',3.1006,'Mwt',18,'MwtH2',2.01588,'Tref',25,'rho',1000,'Vc',2.0681,'Vh',1.9944);
%Cp=specific heat of water, [J/gK];Mwt=mol. wt of H2O, rho=density of
%water/lye[kg/m3],Vc=volume of cold side of heat exchanger[m3],Vh=volume of
%hot side of heat exchanger[m3]
par.Comp = struct('alpha',0.63,'k',1.62,'Tel',25+273,'Pel',3);
par.Storage = struct('VstoH2',965000,'VstoO2',482500,'PoutH2',19,'PoutO2',19,'TstoH2',25+273.15,'TstoO2',25+273.15,'Rg',8.314e-2,'VdispH2',0.4,'VdispO2',0.4);%VstoH2 and VstoO2 are in litres
par.Tw_in = 10;%inlet temperature of the cooling water in lye circulation heat exchanger

%% Parameters for U-I relationship in Ulleberg's model
par.U = struct([]);
par.TherMo = struct([]);
par.EL = struct([]);
for i =1:N
    par.U(i).r1 = 0.000218155;               %ohm m^2
    par.U(i).r2 = -0.000000425;               %ohm m^2 C^-1
    par.U(i).s = 0.1179375;                   %Vs
    par.U(i).t1 = -0.14529;                %A^-1 m^2
    par.U(i).t2 = 11.794;                 %A^-1 m^2 C^-1
    par.U(i).t3 = 395.68;                 %A^-1 m^2 C^-2
    par.U(i).f1 = 120;                   %mA^2 cm^-4
    par.U(i).f2 = 0.98;                  %dimensionless
    
    %% Parameters for the thermal model
    %par.TherMo(i).Ccw = 0.7e3;     %thermal capacity of cooling water, [J s^-1 C^-1]
    par.TherMo(i).CtS = 1.136e3;                       %specific thermal capacity of electrolyzer i.e. Ct/I, [J/A*C]
    %par.TherMo(i).Rt = .167;                           %C W^-1
    par.TherMo(i).Tcwi = 14.5;                          %inlet water temp, [C]
    par.TherMo(i).hc = 5.5;      %convective heat transfer coefficient W/m^2 C
    par.TherMo(i).A_surf = 0.1;  %specific radition area per kA current per cell, [m^2/kA*Ncell]
    
    %% Parameters for Faraday effeciency calculations
    par.EL(i).Utn = 1.482;               %thermoneutral voltage, [V]
    %par.EL(i).Urev = 1.229;              %reversible cell voltage, [V]
    par.EL(i).nc = 230;                   %no. of cells
    par.EL(i).A = 2.6;                   %electrode area of each cell, [m^2]
    par.EL(i).Ta = 20;                   %ambient temp, [C]
    par.EL(i).Tstd = 25;                 %standard temperature, [C]
    par.sigma = 5.672*10^-8; %stefan-boltzmann constant [W/m^2 K^4]
    par.em = 0.8; %emissivity [-]
end
% %El #2, Best performing electrolyzer 
% par.U(2).r1 = par.U(2).r1*0.85;               %ohm m^2
% par.U(2).s = par.U(2).s*0.9;                   %V
% par.U(2).f1 = par.U(2).f1*0.9;                   %mA^2 cm^-4
% par.U(2).f2 = 0.99;
% 
% %El #3, Worst performing electrolyzer
% par.U(3).r1 = par.U(3).r1*1.4;               %ohm m^2
% par.U(3).s = par.U(3).s*1.1;                   %V
% par.U(3).f1 = par.U(3).f1*1.1;                   %mA^2 cm^-4
% par.U(3).f2 = 0.96;

par.N=N;
end