function [lbz,lbx,lbu,ubz,ubx,ubu] = decision_varbound(N)

%constraints on algebriac states
lbu_k = 0*ones(N,1);            %lower bound on cell voltage
ubu_k = inf*ones(N,1);

lbi_k = zeros(N,1);             %lower bound on current
ubi_k = inf*ones(N,1);

lbP_k = 0*ones(N,1);            %lower bound on power
ubP_k = inf*ones(N,1);

lbFeff_k = 0*ones(N,1);         %lower bound on faraday efficiency
ubFeff_k = 1*ones(N,1);

lbnH2_k = zeros(N,1);           %lower bound on hydrogen production
ubnH2_k = inf*ones(N,1);

lbqH2Oloss_k = zeros(N,1);      %lower bound on water loss
ubqH2Oloss_k = inf*ones(N,1);

%constraints on the differential states 
lbT_k = 25*ones(N,1);           %lower bound on the electrolyzer temperature
ubT_k = 80*ones(N,1);

lbMbt = 0*ones(N,1);            %lower bound on the mass in the buffer tank
ubMbt = 2000000*ones(N,1);

lbT_bt_out = 0*ones(N,1);       %lower bound on the temperature of lye leaving the buffer tank
ubT_bt_out = inf*ones(N,1);

lbT_el_in = 0*ones(N,1);        %lower bound on the temperature at electrolyzer inlet
ubT_el_in = inf*ones(N,1);

lbT_cw_out = 0*ones(N,1);       %lower bound on the coolant outlet temperature
ubT_cw_out = inf*ones(N,1);

%constraints on the inputs
lbU_el_k = zeros(N,1);          %lower bound on the electrolyzer voltage
ubU_el_k = inf*ones(N,1);

% lbq_lye_k = 500*ones(N,1);      %lower bound on the lye flowrate
% ubq_lye_k = 10000*ones(N,1);

lbq_lye_k = 6.63e3*ones(N,1);      %lower bound on the lye flowrate
ubq_lye_k = 6.63e3*ones(N,1);

lbq_cw = 1e-2*ones(N,1);        %lower bound on the coolant flow rate
ubq_cw = (80000/3)*ones(N,1);

lbqH2O = 0*ones(N,1);           %lower bound on total water lost during electrolysis
ubqH2O = inf*ones(N,1);

lbz = [lbu_k;lbi_k;lbP_k;lbFeff_k;lbnH2_k;lbqH2Oloss_k];        %lower bounds on all the algebriac variables
lbx = [lbT_k;lbMbt;lbT_bt_out;lbT_el_in;lbT_cw_out];            %lower bounds on all the differential variables
lbu = [lbU_el_k;lbq_lye_k;lbq_cw;lbqH2O];                       %lower bounds on all the input variables
ubz = [ubu_k;ubi_k;ubP_k;ubFeff_k;ubnH2_k;ubqH2Oloss_k];        %upper bounds on all the algebriac variables
ubx = [ubT_k;ubMbt;ubT_bt_out;ubT_el_in;ubT_cw_out];            %upper bounds on all the differential variables
ubu = [ubU_el_k;ubq_lye_k;ubq_cw;ubqH2O];                       %upper bounds on all the input variables

end