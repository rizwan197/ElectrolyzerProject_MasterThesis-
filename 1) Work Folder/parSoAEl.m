function par = parSoAEl(N,T)

par.F = 96485;
par.z = 2;
par.Ncell = 230;
par.MwH2 = 2;
par.A = 2.3; %m2
par.delG = 237*10^3; %J/mol
par.lhvH2 = 3*10^3; %Wh/Nm3
par.polyExp = 1.62; %polytropic exponent for diatomic H2 and O2
par.CompEff = 0.63; %compressor efficiency
par.Tel = 298; %Inlet gas temperature to the compressor
par.Pel = 3; %Inlet gas pressure to the compressor
par.Psto = 19; %Outlet gas pressure form the compressor
par.R = 8.314; %Universal gas constant J/mol K

for k=1:length(T)
par.Urev(k) = 1.5184 - 1.5421e-3*(273+T(k)) + 9.523e-5*(273+T(k))*log((273+T(k))) + ...
    9.84e-8*(273+T(k))^2; %relation for Urev with T from LeRoy eqn. 58
end

for i=1:N
par.U(i).r1 = 8.05*10^-5;
par.U(i).r2 = -2.5*10^-7;
par.U(i).s = 0.185;
par.U(i).t1 = -0.1002;
par.U(i).t2 = 8.424;
par.U(i).t3 = 247.3;

par.U(i).f1 = 250;
par.U(i).f2 = 1;
end

par.U(2).r1 = 8.05*10^-5*0.6;
par.U(2).s = 0.185*0.9;
%par.U(2).t1 = -0.1002*1.1;
par.U(2).f1 = 250*0.9;
par.U(2).f2 = 0.97;

par.U(3).r1 = 8.05*10^-5*1.4;
par.U(3).s = 0.185*1.1;
%par.U(3).t1 = -0.1002*0.9;
par.U(3).f1 = 250*1.1;
par.U(3).f2 = 0.95;

end

