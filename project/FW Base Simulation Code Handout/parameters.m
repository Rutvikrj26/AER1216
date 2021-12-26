% AER1216 Fall 2021 
% Fixed Wing Project Code
%
% parameters.m
%
% Initialization file which generates and stores all required data into the 
% structure P, which is then stored in the workspace. Simulink model calls 
% on this function at the start of every simulation. Code structure adapted
% from Small Unmanned Aircraft: Theory and Practice by R.W. Beard and T. W. 
% McLain. 
% 
% Inputs: 
% N/A
%
% Outputs:
% P                 structure that contains all aerodynamic, geometric, and
%                   initial condition data for the aircraft and simulation.
%
% Last updated: Pravin Wedage 2021-11-09

%% TA NOTE
% An easy way to store parameters for use in simulink is through the use of
% a structure. For example, P.g = 9.81 stores the value of gravitational
% acceleration in the field g that is contained within the structure P.
% Anytime P is called anywhere in the simulation code, the value of P.g is
% accessible. 

%% Parameter Computation
% Initial Conditions
clear all
% compute trim conditions            
P.Va0 = 15;         % initial airspeed (also used as trim airspeed)
P.Va_trim = 28.35; 
P.Va = 28.35;


P.gravity = 9.81;
P.g = 9.81; 

% Aerosonde UAV Data
% physical parameters of airframe

% aerodynamic coefficients

% Control Input limits 
P.delta_e_max = deg2rad(45); % assumed symmetric
P.delta_a_max = deg2rad(45); 
P.delta_r_max = deg2rad(25);

% Initial Conditions % connects with aircraft_dynamics.m, do not modify
% structure field names
P.pn0    = 0;  % initial North position
P.pe0    = 0;  % initial East position
P.pd0    = -2000;  % initial Down position (negative altitude)
P.u0     = P.Va; % initial velocity along body x-axis
P.v0     = 0;  % initial velocity along body y-axis
P.w0     = 0;  % initial velocity along body z-axisu_trim
P.phi0   = 0;  % initial roll angle
P.theta0 = 0;  % initial pitch angle
P.psi0   = 0;  % initial yaw angle
P.p0     = 0;  % initial body frame roll rate
P.q0     = 0;  % initial body frame pitch rate
P.r0     = 0;  % initial body frame yaw rate
P.delta_e0 =0;
P.delta_a0 =0;
P.delta_r0 =0;
P.delta_t0 =0;
                         
P.CLo=0.28;
P.CDo=0.03;
P.Cmo=-0.02338;
P.CLa=3.45;
P.CDa=0.30;
P.Cma=-0.38;
P.CLq=0;
P.CDq=0;
P.Cmq=-3.6;
P.CLdelta_e=-0.36;
P.CDdelta_e=0;
P.Cmdelta_e=-0.5;
P.epsilon=0.1592;

% P.Cxu  = 0;% CTu = -0.096; CDu = 0; Cxu = CTu - CDu
% P.Cxo= -CDo;
% P.Cxa  = -CDa+CLo; 		% C_{x_\alpha}
% P.Cxq  =  -CDq;
% P.Cxdu =  0.0;			% C_{x_{\dot u}}
% P.Cxda =  0.0; 			% C_{x_{\dot \alpha}}
% P.Cxdq =  0.0;

% P.Czu  = 0; %Cz =-CL
% P.Czo = -CLo;
% P.Cza  = -(CDo + CLa);
% P.Czq  = -CLq;
% P.Czdu =  0.0;
% P.Czda =  0;
% P.Czdq =  0.0;

% P.Cmu  =-0.02338  ;
% P.Cma  = -0.38;
% P.Cmq  = -3.6 ;
% P.Cmdu =  0.0;
% P.Cmda = 0;
% P.Cmdq =  0.0;

% P.Cxdelta_e = -CDdelta_e; %C_{x_{\delta_e}}   
% P.Czdelta_e = -CLdelta_e;
% P.Cmdelta_e = -0.5;
% P.Cxdelta_p = 0;  % this and the following delta_p are arbitrary
% P.Czdelta_p = 0;
% P.Cmdelta_p = 0;

P.Cyo=0;
P.Clo=0;
P.Cno=0;

P.Cyb = -0.98;          % C_{y_\beta}
P.Clb = -0.12;
P.Cnb = 0.25;

P.Cyp = 0;
P.Clp = -0.26;
P.Cnp = 0.022;

P.Cyr = 0;
P.Clr =  0.14;
P.Cnr = -0.35;

P.Cydelta_a = 0;          %C_{y_{\delta_a}}
P.Cldelta_a = 0.08;
P.Cndelta_a = 0.06;
P.Cydelta_r = -0.17;
P.Cldelta_r = 0.105;
P.Cndelta_r = -0.032;