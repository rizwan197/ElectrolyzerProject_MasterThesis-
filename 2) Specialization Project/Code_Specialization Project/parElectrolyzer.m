function par = parElectrolyzer(N)

%This script defines values of the input parameters for all electrolyzers.

par.Const = struct('ze',2,'FC',96485,'R',8.314,'Cp',4.18,'CpLye',3.1,'Mwt',18,'MwtH2',2.01588,'Tref',25);%Cp=specific heat of water, [J/gK];Mwt=mol. wt of H2O
par.Comp = struct('alpha',0.63,'k',1.62,'Tel',25+273,'Pel',3);
par.Storage = struct('VstoH2',9500,'VstoO2',4750,'PoutH2',19,'PoutO2',19,'TstoH2',25+273,'TstoO2',25+273,'Rg',8.314e-2,'VdispH2',0.4,'VdispO2',.4);

%% Parameters for U-I relationship in Ulleberg's model
par.U = struct([]);
par.TherMo = struct([]);
par.EL = struct([]);
for i =1:N
    par.U(i).r1 = 8.05e-5;               %ohm m^2
    par.U(i).r2 = -2.5e-7;               %ohm m^2 C^-1
    par.U(i).s = .185;                   %Vs
    par.U(i).t1 = -.1002;                %A^-1 m^2
    par.U(i).t2 = 8.424;                 %A^-1 m^2 C^-1
    par.U(i).t3 = 247.3;                 %A^-1 m^2 C^-2
    par.U(i).f1 = 250;                   %mA^2 cm^-4
    par.U(i).f2 = 0.96;                  %dimensionless
    
    %% Parameters for Ulleberg's thermal model
    par.TherMo(i).Ccw = 0.7e3;                                     %thermal capacity of cooling water, J s^-1 C^-1
    par.TherMo(i).Ct = 625e3;       %625e3                         %thermal capacity of electrolyzer, J C^-1
    par.TherMo(i).tauT = par.TherMo(i).Ccw*par.TherMo(i).Ct;       %thermal time constant, Rt*Ct, [s]
    par.TherMo(i).Rt = .167;                                       %C W^-1
    par.TherMo(i).Tcwi = 14.5;                                     %inlet water temp, [C]
    par.TherMo(i).hcond = 7;                                       %W C^-1
    par.TherMo(i).hconv = 0.02;                                    %W C^-1 per A
    
    %% Parameters for Faraday effeciency calculations
    par.EL(i).Utn = 1.482;               %thermoneutral voltage, [V]
    par.EL(i).Urev = 1.229;              %reversible cell voltage, [V]
    par.EL(i).nc = 21;                   %no. of cells
    par.EL(i).A = .25;                   %electrode area of each cell, [m^2]
    par.EL(i).Ta = 20;                   %ambient temp, [C]
    par.EL(i).Tstd = 25;                 %standard temperature, [C]
end
%El #2
par.U(2).r1 = 8.05e-5*0.85;               %ohm m^2
par.U(2).s = .185*0.9;                   %V
par.U(2).f1 = 250*0.9;                   %mA^2 cm^-4
par.U(2).f2 = 0.97;

%El #3
par.U(3).r1 = 8.05e-5*1.4;               %ohm m^2
par.U(3).s = .185*1.1;                   %V
par.U(3).f1 = 250*1.1;                   %mA^2 cm^-4
par.U(3).f2 = 0.95;

par.N=N;
end