close all;
clear;
clc;

%%%%%%%%%%%%%% Control System %%%%%%%%%%%%%%%%%%
tCycle = 0.001; % [s]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constant Supply Pessure Settings
p_S=180;                    % [bar]     Supply pressure
p_R=1.1;                    % [bar]     Return pressure
%% Initial Settings 
% Initial Values
x_cyl_0 = 0.05;        % [m]
% Cylinder Velocity
v_cyl_0=0;                  % [m/s]     Initial cylinder speed
% System Pressures
p0_B=1.1;                  % [bar]     Initial pressure rod side chamber
p0_A=70.5441;                  % [bar]     Initial pressure piston side chamber
p0_A1=0.99958;                    % [bar]     Initial pressure ocv-dcv volume
p0_C=26.3515;                    % [bar]     Initial pressure pcomp-dcv volume
%% Valve settings
DB = -0.114;         % Valve Deadband
uS_sp = -0.55;       % Valve opening Set Point [-1....1] 
uS0 = -0.00;         % Initial Valve Opening [-1....1] 
%% Volumes
V_CH_A = 1.0e-03;
V_CH_A1 = 0.1e-03;
V_CH_B = 1.25e-03; 
%% Global parameters
k_tanh=0.005;                 % [-]       Constant for steepening tanh function (sign)
%% Over-Centre Valve parameters
p_ocv_set=150;               % [bar]     Crack pressure of the relief valve
mu_ocv=3;                   % [-]       Pilot ratio
k_ocv=15;                   % [l/min/sqrt(bar)] Flow coefficient
p_ocv_fo=350;               % [bar]     Pressure to fully open
tau_ocv=0.0889*0+0.055;             % [s]       Time constant for valve opening dynamics
%% Check Valve parameters (parallel to OCV)
p_cv_set=2.8;               % [bar]     Check valve crash pressure
p_cv_fo=1;                  % [bar]     Additional pressure to fully open
k_cv=15;                    % [l/min/sqrt(bar)] Flow coefficient
%% Pressure Compensator parameters
p_pc_set=7;                 % [bar]     Compensator spring set pressure
tau_pc= 0.009;             % [-]       Time constant pc dynamics
%% Directional Control Valve parameters
% Valve Dynamics
wn_dcv=30;                  % [rad/s]   Natural frequency dcv spool
dr_dcv=0.7;                 % [-]       Damping ratio dcv spool
% Valve Opening edges
OL_pa=-0.114;                 % [-]       Overlap edge p-a
OL_pb=0.114;                  % [-]       Overlap edge p-b
UL_at=-0.05;                % [-]       Underlap edge a-t
UL_bt=0.05;                 % [-]       Underlap edge b-t
% DCV Flow Characteristics
k_DCV_PA= 7.9;           % [(L/min)/(bar^0.5)] (PUMP TO ACTUATOR A-PORT)
k_DCV_PB= 7.78;           % [(L/min)/(bar^0.5)] (PUMP TO ACTUATOR B-PORT)
k_DCV_T= 7.9;                   % [(L/min)/(bar^0.5)]  (ACTUATOR TO TANK)
%% Oil parameters
K_oil = 12000;                % [bar]     Oil Bulk modulus
%% Effective Bulk Modulus ON/OFF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Eff_Bulk_On = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
patm = 1.013e5;                                 % [Pa] Athmospheric Pressure
eta_air = 0.7e-2;                 % [%] Percentage Air in Fluid
kappa_air = 1.4;                                % [%] Adiabatic Constant
beta_F = 12000e5;                               % [Pa] Max. Effective Bulk Modulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cylinder parameters
Bore_diameter=65;           % [mm]      Cylinder piston side diameter
Rod_diameter=35;            % [mm]      Cylinder rod side diameter
A_ratio=1.41;               % [-]       Cylinder area ratio
% Cyl areas
Bore_area=pi/4*Bore_diameter^2;  % [mm^2]
Rod_area=pi/4*Rod_diameter^2;    % [mm^2]
Annulus_area=Bore_area-Rod_area; % [mm^2]
% Common Parameters
Stroke=500;                      %[mm]   stroke
cyl_KLi=0;                       %[m^3/s/Pa] coefficient of internal leakage
cyl_x_compen=0;                  %[m]    origin compensation
% Other Parameters (from validation SCC)
cyl_c_visc=1500;                 % [kg/s]  Coeff. viscous friction (15000)
cyl_F_coul=75;                    % [N]     Coulomb force (75)
cyl_FH_stat=500;                 % [N]     Static friction force (1450)
cyl_CH_stat=1/50;                 % [s/m]   Static friction time constant  
cyl_gamma=250;                    % [-]     Approximation for friction force
% Cylinder End stop
cyl_k_end=6e9;                   %[N/m]  spring constant used for stroke limits on cylinder
cyl_c_end=9e8;  


