%This script defines values for all the input parameters.

parU = struct('r1',8.05e-5,'r2',-2.5e-7,'s',.185,'t1',-.1002,'t2',8.424,'t3',247.3,'f1',250,'f2',0.96);
parTherMo = struct('tauT',29*3600,'Ccw',0.7e3,'Ct',625e3,'Rt',0.167,'Tcwi',14.5,'hcond',7,'hconv',0.02);
parEL = struct('Utn',1.482,'Urev',1.229,'nc',21,'A',.25,'Ta',20,'Tstd',25);
parConst = struct('ze',2,'FC',96485,'R',8.314);
parComp = struct('alpha',0.63,'k',1.62,'Tel',25+273,'Pel',3);
parStorage = struct('VstoH2',5,'VstoO2',5,'PoutH2',19,'PoutO2',19,'TstoH2',25+273,'TstoO2',25+273,'Rg',8.314e-2,'VdispH2',.4,'VdispO2',.4);

%% Parameters for U-I relationship in Ulleberg's model
r1 = parU.r1;               %ohm m^2
r2 = parU.r2;               %ohm m^2 C^-1
s = parU.s;                 %V
t1 = parU.t1;               %A^-1 m^2
t2 = parU.t2;               %A^-1 m^2 C^-1
t3 = parU.t3;               %A^-1 m^2 C^-2
f1 = parU.f1;               %mA^2 cm^-4
f2 = parU.f2;               %dimensionless

%% Parameters for Ulleberg's thermal model
tauT = parTherMo.tauT;      %thermal time constant, Rt*Ct, [s]
Ccw = parTherMo.Ccw;        %J s^-1 C^-1
Ct = parTherMo.Ct;          %J C^-1
Rt = parTherMo.Rt;          %C W^-1
Tcwi = parTherMo.Tcwi;      %inlet water temp, [C]
hcond = parTherMo.hcond;    %W C^-1
hconv = parTherMo.hconv;    %W C^-1 per A

%% Parameters for Faraday effeciency calculations
Utn = parEL.Utn;            %thermoneutral voltage, [V]
Urev = parEL.Urev;          %reversible cell voltage, [V]
nc = parEL.nc;              %no. of cells
A = parEL.A;                %electrode area of each cell, [m^2]
Ta = parEL.Ta;              %ambient temp, [C]
Tstd= parEL.Tstd;           %standard temperature, [C]

%% Constants
ze = parConst.ze;           %number of electrons transferred per reaction 
FC = parConst.FC;           %faraday constant 
R = parConst.R;             %gas constant, [J mol^-1 K^-1,]

%% Parameters for compressor calculations
alpha = parComp.alpha;      %compressor efficient, 63%
k = parComp.k;              %polytropic exponent, 1.62
Tel = parComp.Tel;          %temperature of inlet to compressor after cooler, [Kelvin]
Pel = parComp.Pel;          %pressure of electrolyzer outlet, [bar]

%% Parameters for storage system model
VstoH2 = parStorage.VstoH2;     %volume of H2 storage, litres
VstoO2 = parStorage.VstoO2;     %volume of O2 storage,litres
PoutH2 = parStorage.PoutH2;     %outlet H2 pressure, bar
PoutO2 = parStorage.PoutO2;     %outlet O2 pressure, bar
TstoH2 = parStorage.TstoH2;     %temp of H2 storage, kelvin
TstoO2 = parStorage.TstoO2;     %temp of H2 storage, kelvin
Rg = parStorage.Rg;             %gas constant in l bar K^-1 mol^-1
VdispH2 = parStorage.VdispH2;   %steady state valve displacement for hydrogen outlet
VdispO2 = parStorage.VdispO2;   %steady state valve displacement for oxygen outlet