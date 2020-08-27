# MasterThesis
Branches of this github repository corresponds to the different flowsheet designs. Key is defined as below:
1.	master - This file includes the work in progress on control structure design

2. CtrlStrDesign_v2SOC_Correct - This includes the work on the design of the supervisory control layer and is ongoing work on the selection of SOC 

3. SystemAnalysis - Code of electrolyzer plant flowsheet with shared BoP and power systems. The electrolyzer stacks are arranged electrically in parallel. This branch consists of the following files: 
3.a) "main.m" - This is the main code file, run it for steady state optimization. However, put a break at line no. 122 to stop the dynamic simulation when only steady state optimization problem is solved. 
3.b) "parElectrolyzer.m" - Use this file to change the heat exchanger size, state of the electrolyzers (like new or old) for a given flowsheet. 
3.c) "model.m" - This file includes the dynamic model for the electrolyzer plant flowsheet. Casadi framework is used for symbolic modelling in MATLAB. 
3.d) "El_SteadyStateOptimization" - This file solves the steady state optimization problem. All the ODE equations for dyanmic states in the model are equated to 0. IPOPT solver from SUNDIALS suite is used for solving the optimization problem. The bounds on the all the decision variables for the optimization problem are also included in this file.

4.	Separated_Electrolyzers - Code of electrolyzer plant with separate BoP and power systems. This branch consists of the following files: 
4.a) "main.m" - This is the main code file, run it for steady state optimization. However, put a break at line no. 102 to stop the dynamic simulation when only steady state optimization problem is solved. 
4.b) "parElectrolyzer.m" - Use this file to change the heat exchanger size, state of the electrolyzers (like new or old) for a given flowsheet. 
4.c) "model.m" - This file includes the dynamic model for the electrolyzer plant flowsheet. Casadi framework is used for symbolic modelling in MATLAB. 
4.d) "El_SteadyStateOptimization" - This file solves the steady state optimization problem. All the ODE equations for dyanmic states in the model are equated to 0. IPOPT solver from SUNDIALS suite is used for solving the optimization problem. 
4.e) "decision_varbound.m" - This file includes the bounds that are defined on all the decision variables for optimization problem (i.e. dynamic state variables, algebriac variables and inputs variables)

