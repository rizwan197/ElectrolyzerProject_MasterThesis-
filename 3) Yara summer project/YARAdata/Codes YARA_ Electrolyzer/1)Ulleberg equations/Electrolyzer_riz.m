clc
clear all
close all

parU = struct('r1',[8.05]*10^-5,'r2',-[2.5]*10^-7,'s',[0.185],'t1',-[0.1002],'t2',[8.424],'t3',[247.3]); %I-U curve parameters
F = 96485;
z = 2;
Ncell = 21;
delG = 237*10^3; %kJ/mol, at standard conditions, 25C and 1 bar
A = .25; %m2
I = [0:0.1:875];%in A (I/A=350, i.e. I=350*.25/.1)
T = [25]; %in Celsius

for i=1:length(T)
    Urev(i) = 1.5184 - 1.5421e-3*(273+T(i)) + 9.523e-5*(273+T(i))*log((273+T(i))) + 9.84e-8*(273+T(i))^2; %relation for Urev with T from LeRoy eqn. 58
end

%nH2 vs current, Fig. 9, Ulleberg's paper
parF = struct('f1',[250],'f2',[0.9 0.96 1]);

u=[];
% I-U curve, Fig. 5, Ulleberg's paper
for i=1:length(T)
    for j=1:length(parU.r1)
        for k=1:length(parU.r2)
            for l=1:length(parU.s)
                for m=1:length(parU.t1)
                    for n=1:length(parU.t2)
                        for o=1:length(parU.t3)
                            for p=1:length(parF.f1)
                                for q=1:length(parF.f2)
                                    U = Urev(i) + (((parU.r1(j)+parU.r2(k)*T(i)).*I)./A) + parU.s(l)*log10(((parU.t1(m)+(parU.t2(n)/T(i))+(parU.t3(o)/T(i)^2)).*I/A)+1);
                                    u=[u U];
                                    figure(1)
                                    plot(0.1*I/A,U)% I-U curve of an electrolyser is in mA/cm2 vs V
                                    xlabel('Current Density, mA/cm^2')
                                    ylabel('Voltage, V/cell')
                                    legend(strcat('T=',num2str(T'),', r_1=',num2str(parU.r1'),', r_2=',num2str(parU.r2'),', s=',num2str(parU.s'),', t_1=',num2str(parU.t1'),', t_2=',num2str(parU.t2'),', t_3=',num2str(parU.t3')),'location','southeast')
                                    hold on
                                    
                                    %nH2 vs current, Fig. 9, Ulleberg's paper
                                    Feff = (((0.1*I./A).^2)./(parF.f1(p)+(0.1*I./A).^2)).*parF.f2(q); %to match the dimensions for f1 units, I/A is multiplied to 0.1
                                    rateH2 = Feff.*Ncell.*I./(z*F); %in mol/s
                                    vstd = 0.0224136;
                                    nH2 = rateH2*vstd*3600;
                                    figure(2)
                                    plot(I,nH2)
                                    xlabel('Current, A')
                                    ylabel('H_2 Flow Rate, Nm^3/hr')
                                    legend(strcat('f_1=',num2str(parF.f1'),', f_2=',num2str(parF.f2')),'location','southeast')
                                    hold on
                                    
                                    %plot between efficiency percentage and current density, Fig. 6, Ulleberg's paper
                                    figure(3)
                                    plot(0.1*I/A,100*Feff)
                                    xlabel('Current Density, mA/cm^2')
                                    ylabel('Efficiency, %')
                                    legend(strcat('f_1=',num2str(parF.f1'),', f_2=',num2str(parF.f2')),'location','southeast')
                                    hold on
                                    
                                    %plot between specific electricity and rate of H2 production, Vidar's plot
                                    figure(4)
                                    E=(U.*I)*3600;%in Wh
                                    SpElec = 1e-6*E./rateH2; %MWh/(mol/hr)
                                    plot(rateH2*3600,SpElec)
                                    ylim([40 80])
                                    ylabel('Specific electricity, MWh/mol H_2')
                                    xlabel('Molar flow rate of H_2, mol/hr')
                                    hold on
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
legend(strcat('T=',num2str(T'),', r_1=',num2str(parU.r1'),', r_2=',num2str(parU.r2'),', s=',num2str(parU.s'),', t_1=',num2str(parU.t1'),', t_2=',num2str(parU.t2'),', t_3=',num2str(parU.t3'),', f_1=',num2str(parF.f1'),', f_2=',num2str(parF.f2')),'location','southeast')

