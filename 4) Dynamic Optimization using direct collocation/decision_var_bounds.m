function[lbx,lbz,lbu,ubx,ubz,ubu] = decision_var_bounds(N)

par = parElectrolyzer(N);
%constraints on states
lbu_k = 0*ones(par.N,1);%lower bound on cell voltage
ubu_k = inf*ones(par.N,1);

lbi_k = zeros(par.N,1);%lower bound on current
ubi_k = inf*ones(par.N,1);

lbP_k = 0*ones(par.N,1);%lower bound on power
ubP_k = inf*ones(par.N,1);

lbFeff_k = 0*ones(par.N,1);%lower bound on faraday efficiency
ubFeff_k = 1*ones(par.N,1);

lbnH2_k = zeros(par.N,1);%lower bound on hydrogen production
ubnH2_k = inf*ones(par.N,1);

lbqH2Oloss_k = zeros(par.N,1);%lower bound on water loss
ubqH2Oloss_k = inf*ones(par.N,1);

lbnH2el_net = 0;%lower bound on net hydrogen production from the electrolyzer
ubnH2el_net = inf;

lbnH2out_net = 0;%lower bound on hydrogen from the storage tank outlet 
ubnH2out_net = inf;

lbnO2el_net = 0;%lower bound on net oxygen production from the electrolyzer
ubnO2el_net = inf;

lbnO2out_net = 0;%lower bound on oxygen from the storage tank outlet
ubnO2out_net = inf;

lbT_el_out = 0;%lower bound on the temperature at electrolyzer outlet
ubT_el_out = inf;

lbT_k = 25*ones(par.N,1);%lower bound on the electrolyzer temperature
ubT_k = 80*ones(par.N,1);

lbPstoH2 = 20;%lower bound on the hydrogen storage pressure
ubPstoH2 = 30;

lbPstoO2 = 20;%lower bound on the oxygen storage pressure
ubPstoO2 = 30;

lbMbt = 0;%lower bound on the mass in the buffer tank 
ubMbt = 6000000;

lbT_bt_out = 0;%lower bound on the temperature of lye leaving the buffer tank 
ubT_bt_out = inf;

lbT_el_in = 0;%lower bound on the temperature at electrolyzer inlet
ubT_el_in = inf;

lbT_cw_out = 0;%lower bound on the coolant outlet temperature
ubT_cw_out = inf;

%constraints on the inputs
lbU_el_k = zeros(par.N,1);%lower bound on the electrolyzer voltage
ubU_el_k = inf*ones(par.N,1);
lbq_lye_k = 500*ones(par.N,1);%lower bound on the lye flowrate
ubq_lye_k = 10000*ones(par.N,1);
lbq_cw = 0.01000;%lower bound on the coolant flow rate
ubq_cw = 80000;
lbzH2 = 0;%lower bound on hydrogen outlet valve opening
ubzH2 = 1;
lbzO2 = 0;%lower bound on oxygen outlet valve opening
ubzO2 = 1;
lbqH2O = 0;%lower bound on total water lost during electrolysis
ubqH2O = inf;

lbz = [lbu_k;lbi_k;lbP_k;lbFeff_k;lbnH2_k;lbqH2Oloss_k;lbnH2el_net;...
    lbnH2out_net;lbnO2el_net;lbnO2out_net;lbT_el_out];
lbx = [lbT_k;lbPstoH2;lbPstoO2;lbMbt;lbT_bt_out;lbT_el_in;lbT_cw_out];
lbu = [lbU_el_k;lbq_lye_k;lbq_cw;lbzH2;lbzO2;lbqH2O];

ubz = [ubu_k;ubi_k;ubP_k;ubFeff_k;ubnH2_k;ubqH2Oloss_k;ubnH2el_net;...
    ubnH2out_net;ubnO2el_net;ubnO2out_net;ubT_el_out];
ubx = [ubT_k;ubPstoH2;ubPstoO2;ubMbt;ubT_bt_out;ubT_el_in;ubT_cw_out];
ubu = [ubU_el_k;ubq_lye_k;ubq_cw;ubzH2;ubzO2;ubqH2O];
end