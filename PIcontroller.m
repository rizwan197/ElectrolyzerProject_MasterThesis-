function u = PIcontroller(u0,Kc,tauI,e,eint)
%This is a function file for the PI controller
%u0 = initial value of MV
%e = error term, i.e. y - yset
%eint = integrated error for integral action
%Kc = proportional gain
%tauI = integral time constant
u = u0 - Kc*e - (Kc/tauI)*eint;
end