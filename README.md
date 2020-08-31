# MasterThesis
Branches of this github repository corresponds to the different flowsheet designs. Key is defined as below:
1.	master - This includes the work on the design of the supervisory control layer and is ongoing work on the selection of SOC. Don't touch this branch!! 

2. SystemAnalysis - Code of electrolyzer plant flowsheet with shared BoP and power systems. The electrolyzer stacks are arranged electrically in parallel. This branch consists of the following files: 
2.a) "main.m" - This is the main code file, run it for steady state optimization. 
2.b) "parElectrolyzer.m" - Use this file to change the heat exchanger size, state of the electrolyzers (like new or old) for a given flowsheet. 
2.c) "model.m" - This file includes the dynamic model for the electrolyzer plant flowsheet. Casadi framework is used for symbolic modelling in MATLAB. 
2.d) "El_SteadyStateOptimization" - This file solves the steady state optimization problem. All the ODE equations for dyanmic states in the model are equated to 0. IPOPT solver from SUNDIALS suite is used for solving the optimization problem. The bounds on the all the decision variables for the optimization problem are also included in this file.
2.e) "UIcurves.m" - This file includes the matlab code to plot the characteristic curves for the electrolyzers. These plots are included in Chapter 3 (Fig 3.5) of the master thesis. 
2.f) Folder "Active region map plotfiles" - This folder includes the steady state optimization results (scaled to limiting values, so that all variable values are between 0-1) for all the different flowsheet designs. This data is then used to plot active regions in different power ranges, i.e. low, medium and high. The plots genreated by plotfile "ACplotfile.m" in this folder are included in Chapter 5 (Fig 5.5) of the master thesis.
2.g) Folder "State1 El data files" and "State2 El data files" - These folders include data files for all the flowsheet's steady state optimization results. The variable values are not scaled and loss formula is used to compare different flowsheet designs w.r.t most capital intensive flowsheet design. MATLAB script "FScompPlotfileS1.m" and "FScompPlotfileS2.m" generate the figures in Chapter 4 (Fig. 4.5, 4.6, 4.7, 4.8 & 4.9) of the master thesis.

3.	Separated_Electrolyzers - Code of electrolyzer plant with separate BoP and power systems. This branch consists of the following files: 
3.a) "main.m" - This is the main code file, run it for steady state optimization. 
3.b) "parElectrolyzer.m" - Use this file to change the heat exchanger size, state of the electrolyzers (like new or old) for a given flowsheet. 
3.c) "model.m" - This file includes the dynamic model for the electrolyzer plant flowsheet. Casadi framework is used for symbolic modelling in MATLAB. 
3.d) "El_SteadyStateOptimization" - This file solves the steady state optimization problem. All the ODE equations for dyanmic states in the model are equated to 0. IPOPT solver from SUNDIALS suite is used for solving the optimization problem. 
3.e) "decision_varbound.m" - This file includes the bounds that are defined on all the decision variables for optimization problem (i.e. dynamic state variables, algebriac variables and inputs variables)

4. RegulatoryLayerControl - This branch includes the code for the implemetation of the regualtory control layer to the flowsheet with shared BoP systems and heat exchanger design HX_EoL (i.e. based on the cooling requirement of the electrolyzer plant at the end of the lifetime). There are following files in this branch:
4.a) "main.m" - This is main code file, all the controllers in the regulatory layer are tuned and the plant is simulated for 3 hours. The disturbance is a step in q_cw (i.e. cooling water flowrate), controllers are able to maintain the states at the desired setpoint values.
4.b) "modelPI.m" - For implementaion of the PI controllers, additional states are created for intergrated error term (i.e. eint) therefore the old model file is updated and renamed as "modelPI.m". This model file is used for dyanmic simulation when the controllers are implemented.
4.c) "mainOld.m" - This is the old main file, this doesn't include any controllers. This file is to be used along with model.m for open loop simulation to the step changes for controller tuning.  
4.d) "model.m" - This file includes the dynamic model for the electrolyzer plant flowsheet. Casadi framework is used for symbolic modelling in MATLAB. The additional states for integrated error term aren't present here in this model file. Also, this model file can be used directly if we choose to implement all the PID contoller in the velocity form (which doesn't need these additional states for integrated error terms)
4.e) "parElectrolyzer.m" - Use this file to change the heat exchanger size, state of the electrolyzers (like new or old) for a given flowsheet. 
4.f) "PIcontroller.m" - Function file for PI controller equation.
4.g) "El_SteadyStateOptimization" - This file solves the steady state optimization problem. All the ODE equations for dyanmic states in the model are equated to 0. IPOPT solver from SUNDIALS suite is used for solving the optimization problem. Use this only when performing step tests for controller tuning to give good initial values for the dynamic simulation.
4.h) Folder "Example PID implementation" - This folder consist of files for finite PID controller implementation. The CSTR example from Sigurd's book is used. The code files for different ways of implementing a discrete PID controller (i.e. with eint term or the velocity form) are included.
