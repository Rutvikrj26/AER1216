function aircraft_dynamics1216_PROJECT(block)
%MSFUNTMPL_BASIC A Template for a Level-2 MATLAB S-Function
%   The MATLAB S-function is written as a MATLAB function with the
%   same name as the S-function. Replace 'msfuntmpl_basic' with the 
%   name of your S-function.
%
%   It should be noted that the MATLAB S-function is very similar
%   to Level-2 C-Mex S-functions. You should be able to get more
%   information for each of the block methods by referring to the
%   documentation for C-Mex S-functions.
%
%   Copyright 2003-2010 The MathWorks, Inc.

% AER1216 Fall 2021 
% Fixed Wing Project Code
%
% aircraft_dynamics.m
%
% Fixed wing simulation model file, based on the Aerosonde UAV, with code
% structure adapted from Small Unmanned Aircraft: Theory and Practice by 
% R.W. Beard and T. W. McLain. 
% 
% Inputs: 
% delta_e           elevator deflection [deg]
% delta_a           aileron deflection [deg]
% delta_r           rudder deflection [deg]
% delta_t           normalized thrust []
%
% Outputs:
% pn                inertial frame x (north) position [m]
% pe                inertial frame y (east) position [m]
% pd                inertial frame z (down) position [m]
% u                 body frame x velocity [m/s]
% v                 body frame y velocity [m/s]
% w                 body frame z velocity [m/s]
% phi               roll angle [rad]
% theta             pitch angle [rad]
% psi               yaw angle [rad]
% p                 roll rate [rad/s]
% q                 pitch rate [rad/s]
% r                 yaw rate [rad/s]
%
% Last updated: Pravin Wedage 2021-11-09

%% TA NOTE
% The code segements you must modify are located in the derivatives
% function in this .m file. Modify other sections at your own risk. 


%
% The setup method is used to set up the basic attributes of the
% S-function such as ports, parameters, etc. Do not add any other
% calls to the main body of the function.
%
setup(block);

end 


%% Function: setup ===================================================
% Abstract:
%   Set up the basic characteristics of the S-function block such as:
%   - Input ports
%   - Output ports
%   - Dialog parameters
%   - Options
%
%   Required         : Yes
%   C-Mex counterpart: mdlInitializeSizes
%
function setup(block)

% Register number of ports
block.NumInputPorts  = 1;
block.NumOutputPorts = 1;

% Setup port properties to be inherited or dynamic
block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;

% Override input port properties
for i = 1:block.NumInputPorts
    block.InputPort(i).Dimensions        = 4;
    block.InputPort(i).DatatypeID  = 0;  % double
    block.InputPort(i).Complexity  = 'Real';
    block.InputPort(i).DirectFeedthrough = false; % important to be false 
end

% Override output port properties
for i = 1:block.NumOutputPorts
    block.OutputPort(i).Dimensions       = 12;
    block.OutputPort(i).DatatypeID  = 0; % double
    block.OutputPort(i).Complexity  = 'Real';
%     block.OutputPort(i).SamplingMode = 'Sample';
end

% Register parameters
block.NumDialogPrms     = 1;
P = block.DialogPrm(1).Data; % must duplicate this line in each function

% Register sample times
%  [0 offset]            : Continuous sample time
%  [positive_num offset] : Discrete sample time
%
%  [-1, 0]               : Inherited sample time
%  [-2, 0]               : Variable sample time
block.SampleTimes = [0 0];

% Register multiple instances allowable
% block.SupportMultipleExecInstances = true;

% Register number of continuous states
block.NumContStates = 12;

% Specify the block simStateCompliance. The allowed values are:
%    'UnknownSimState', < The default setting; warn and assume DefaultSimState
%    'DefaultSimState', < Same sim state as a built-in block
%    'HasNoSimState',   < No sim state
%    'CustomSimState',  < Has GetSimState and SetSimState methods
%    'DisallowSimState' < Error out when saving or restoring the model sim state
block.SimStateCompliance = 'DefaultSimState';

% -----------------------------------------------------------------
% The MATLAB S-function uses an internal registry for all
% block methods. You should register all relevant methods
% (optional and required) as illustrated below. You may choose
% any suitable name for the methods and implement these methods
% as local functions within the same file. See comments
% provided for each function for more information.
% -----------------------------------------------------------------

% block.RegBlockMethod('PostPropagationSetup',    @DoPostPropSetup); % discrete states only
block.RegBlockMethod('SetInputPortSamplingMode', @SetInpPortFrameData);
block.RegBlockMethod('InitializeConditions',    @InitializeConditions);
% block.RegBlockMethod('Start',                   @Start); % Initialize Conditions is used
block.RegBlockMethod('Outputs',                 @Outputs); % Required
% block.RegBlockMethod('Update',                  @Update); % only required for discrete states
block.RegBlockMethod('Derivatives',             @Derivatives); % Required for continuous states
block.RegBlockMethod('Terminate',               @Terminate); % Required

end 


%% PostPropagationSetup:
%   Functionality    : Setup work areas and state variables. Can
%                      also register run-time methods here
%   Required         : No
%   C-Mex counterpart: mdlSetWorkWidths
%
function DoPostPropSetup(block)
block.NumDworks = 1;
  
  block.Dwork(1).Name            = 'x1';
  block.Dwork(1).Dimensions      = 1;
  block.Dwork(1).DatatypeID      = 0;      % double
  block.Dwork(1).Complexity      = 'Real'; % real
  block.Dwork(1).UsedAsDiscState = true;

end


%% InitializeConditions:
%   Functionality    : Called at the start of simulation and if it is 
%                      present in an enabled subsystem configured to reset 
%                      states, it will be called when the enabled subsystem
%                      restarts execution to reset the states.
%   Required         : No
%   C-MEX counterpart: mdlInitializeConditions
%
function InitializeConditions(block)

% Rename parameters
P = block.DialogPrm(1).Data; % must duplicate this line in each function

% Initialize continuous states
block.ContStates.Data(1) = P.pn0; 
block.ContStates.Data(2) = P.pe0;
block.ContStates.Data(3) = P.pd0;
block.ContStates.Data(4) = P.u0;
block.ContStates.Data(5) = P.v0;
block.ContStates.Data(6) = P.w0;
block.ContStates.Data(7) = P.phi0;
block.ContStates.Data(8) = P.theta0;
block.ContStates.Data(9) = P.psi0;
block.ContStates.Data(10) = P.p0;
block.ContStates.Data(11) = P.q0;
block.ContStates.Data(12) = P.r0;

end 

%% Start:
%   Functionality    : Called once at start of model execution. If you
%                      have states that should be initialized once, this 
%                      is the place to do it.
%   Required         : No
%   C-MEX counterpart: mdlStart
%
function Start(block)

block.Dwork(1).Data = 0;

end 

%% Input Port Sampling Method:
function SetInpPortFrameData(block, idx, fd)
  
  block.InputPort(idx).SamplingMode = 'Sample';
  for i = 1:block.NumOutputPorts
    block.OutputPort(i).SamplingMode  = 'Sample';   
  end
end

%% Outputs:
%   Functionality    : Called to generate block outputs in
%                      simulation step
%   Required         : Yes
%   C-MEX counterpart: mdlOutputs
%
function Outputs(block)

temp_mat = zeros(block.NumContStates,1); % thirteen states
for i = 1:block.NumContStates
     temp_mat(i) = block.ContStates.Data(i);
end

block.OutputPort(1).Data = temp_mat; % states

% for i = 1:block.NumOutputPorts
%     block.OutputPort(1).Data(i) = block.ContStates.Data(i);
% end

end 


%% Update:
%   Functionality    : Called to update discrete states
%                      during simulation step
%   Required         : No
%   C-MEX counterpart: mdlUpdate
%
function Update(block)

block.Dwork(1).Data = block.InputPort(1).Data;

end 


%% Derivatives:
%   Functionality    : Called to update derivatives of
%                      continuous states during simulation step
%   Required         : No
%   C-MEX counterpart: mdlDerivatives
%
function Derivatives(block)

% Rename parameters
P = block.DialogPrm(1).Data; % must duplicate this line in each function

% compute inertial constants
% K = ;
% k1 = ;
% k2 = ;
% k3 = ;
% k4 = ;
% k5 = ;
% k6 = ;
% k7 = ;
% k8 = ;

% map states and inputs
pn    = block.ContStates.Data(1);
pe    = block.ContStates.Data(2);
h    = block.ContStates.Data(3);
u     = block.ContStates.Data(4);
v     = block.ContStates.Data(5);
w     = block.ContStates.Data(6);
phi   = block.ContStates.Data(7);
theta = block.ContStates.Data(8);
psi   = block.ContStates.Data(9);
p     = block.ContStates.Data(10);
q     = block.ContStates.Data(11);
r     = block.ContStates.Data(12);
delta_e = block.InputPort(1).Data(1)*pi/180 ; % converted inputs to radians
delta_a = block.InputPort(1).Data(2)*pi/180 ; % converted inputs to radians
delta_r = block.InputPort(1).Data(3)*pi/180 ; % converted inputs to radians
delta_t = block.InputPort(1).Data(4);


%%  Basic Inputs

g = 9.81; 			% gravity (m/s^2)
m = 13.5; 	    	% mass (kg)
W = m*g;   % weight (N)
Ixx = 0.8244;   % moment of inertia (kg m^2)
Iyy = 1.135;   % 
Izz = 1.759;   %
Ixz = 0.1204;
Izx = Ixz;

S = 0.55; 			% wing plan area (m^2)
c_bar = 0.18994; 	    % mean chord (m) 
b = 2.8956; 			% wing span (m) 
AR = (b^2)/S;
e = 0.9;     % oswald factor, CHECK
k = 1 /(pi*e*AR);

C_T=0.78; % We need to redefine


%  (2) flight conditions reference point ****

Va_trim=P.Va_trim;
Va = P.Va;
Ue = 28.35        ; 	% steady flight velocity (m/s): 
alpha = 0;
theta_e = 0;
Sprop = 0.2027;

rho = 1.0582; 	% density (kg/m^3) %% CHECK
% =========  dimensional derivatives ================
% q = 0.5*rho*Ue^2;
qs = rho*Ue*S/2;
qs1 = rho*Va*S/2;
qs1_2 = rho*Va*Va*S/2;

%% Longitudinal terms

Cxo = -P.CDo*cos(alpha) + P.CLo*sin(alpha);
Cxa = -P.CDa*cos(alpha) + P.CLa*sin(alpha);
Cxq = -P.CDq*cos(alpha) + P.CLq*sin(alpha);
Cxdelta_e = -P.CDdelta_e*cos(alpha) + P.CLdelta_e*sin(alpha);

Xu  = (2*qs*Cxo - rho*Sprop*Ue*C_T)/m;%
Xw  = Cxa*qs/m;
Xq  = Cxq*qs1*c_bar/(2*m);
Xdelta_e = Cxdelta_e*qs1_2/m;
Xdelta_t = (rho*Sprop*C_T*k^2)/m;

Czo = -P.CLo*cos(alpha) - P.CDo*sin(alpha);
Cza = -P.CLa*cos(alpha) - P.CDa*sin(alpha);
Czq = -P.CLq*cos(alpha) - P.CDq*sin(alpha);
Czdelta_e = -P.CLdelta_e*cos(alpha) - P.CDdelta_e*sin(alpha);

Zu  = 2*qs*Czo/m ;
Zw  = Cza*qs/m;
Zq  = Ue + Czq*qs1*c_bar/(2*m); %% no m term in ss matrix
Zdelta_e = Czdelta_e*qs1_2/m;

Cmo = P.Cmo;
Cma = P.Cma;
Cmq = P.Cmq;
Cmdelta_e = P.Cmdelta_e;

Mu  = 2*qs*c_bar*Cmo/Iyy;
Mw  = Cma * qs * c_bar/Iyy;
Mq  = Cmq*qs1*c_bar*c_bar/(2*Iyy);
Mdelta_e = Cmdelta_e*qs1_2*c_bar/Iyy;

%%  Lateral terms

Cyo = P.Cyo;
Cyb = P.Cyb;
Cyp = P.Cyp;
Cyr = P.Cyr;
Cydelta_a = P.Cydelta_a;
Cydelta_r = P.Cydelta_r;

Yv = 2*Cyo*qs/m + qs*Cyb/m;
Yp = Cyp*qs1*b/(2*m);
Yr = -Ue + (Cyr*qs1*b/(2*m));
Ydelta_a = Cydelta_a*qs1_2/m;
Ydelta_r = Cydelta_r*qs1_2/m;

Cpb = P.Clb;
Cpp = P.Clp;
Cpr = P.Clr;
Cpdelta_a = P.Cldelta_a;
Cpdelta_r = P.Cldelta_r;

Lv = Cpb*qs*b;
Lp = Cpp*qs1*b^2/2;
Lr = Cpr*qs*b^2/2;
Ldelta_a = Cpdelta_a * qs1_2 * b;
Ldelta_r = Cpdelta_r * qs1_2 * b;

Crb = P.Cnb;
Crp = P.Cnp;
Crr = P.Cnr;
Crdelta_a = P.Cndelta_a;
Crdelta_r = P.Cndelta_r;

Nv = Crb * qs * b;
Np = Crp * qs1 * b^2 / 2;
Nr = Crr * qs1 * b^2 / 2;
Ndelta_a = Crdelta_a * qs1_2 * b;
Ndelta_r = Crdelta_r * qs1_2 * b;

%% Final Matrices

X_long=[u-Va w q theta h]';
U_long=[delta_e delta_t]';

Along=[Xu Xw Xq -g*cos(theta_e) 0;
       Zu Zw Zq -g*sin(theta_e) 0;
       Mu Mw Mq 0 0
       0 0 1 0 0;
       sin(theta_e) -cos(theta_e) 0 Ue*cos(theta_e) 0];

Blong=[Xdelta_e Xdelta_t;
       Zdelta_e 0;
       Mdelta_e 0;
       0 0;
       0 0];
        
X_long_dot=Along*X_long + Blong*U_long;

udot=X_long_dot(1);
wdot=X_long_dot(2);
qdot=X_long_dot(3);
thetadot=X_long_dot(4);
hdot = X_long_dot(5);

%Lateral 
X_latr=[v p r phi psi]';
U_latr=[delta_a delta_r]';

Alatr=[Yv Yp Yr g*cos(theta_e) 0; 
       Lv Lp Lr    0           0;
       Nv Np Nr    0           0;
       0 1 tan(theta_e) 0      0;
       0 0 sec(theta_e) 0      0];

Blatr=[Ydelta_a Ydelta_r;
       Ldelta_a Ldelta_r;
       Ndelta_a Ndelta_r;
       0        0       ;
       0        0       ];

Blatr(5,:)=[0 0];

X_latr_dot=Alatr*X_latr + Blatr*U_latr;

vdot=X_latr_dot(1);
pdot=X_latr_dot(2);
rdot=X_latr_dot(3);
phidot=X_latr_dot(4);
psidot=X_latr_dot(5);

C_BE=[cos(theta)*cos(psi) sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi) cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi);
    cos(theta)*sin(psi) sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(psi) cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi);
    -sin(theta) sin(phi)*cos(theta) cos(phi)*cos(theta)]';

X_e=C_BE*[u v w]';
pndot = X_e(1); % X Value
pedot =X_e(2) ; % Y Value

% rotation matrix

% Aerodynamic Coefficients 
% compute the nondimensional aerodynamic coefficients here

% aerodynamic forces and moments
% compute the aerodynamic forces and moments here

% propulsion forces and moments
% compute the propulsion forces and moments here

% gravity
% compute the gravitational forces here

% total forces and moments (body frame)

% state derivatives
% the full aircraft dynamics model is computed here
% pdot = ;
% pndot = ;
% pedot = ;
% pddot = ;

% udot = ;
% vdot = ;
% wdot = ;
% 
% phidot = ;
% thetadot = ;
% psidot = ;
% 
% pdot = ;
% qdot = ;
% rdot = ;

% map derivatives
block.Derivatives.Data(1) = pndot;
block.Derivatives.Data(2) = pedot;
block.Derivatives.Data(3) = hdot;
block.Derivatives.Data(4) = udot;
block.Derivatives.Data(5) = vdot;
block.Derivatives.Data(6) = wdot;
block.Derivatives.Data(7) = phidot;
block.Derivatives.Data(8) = thetadot;
block.Derivatives.Data(9) = psidot;
block.Derivatives.Data(10)= pdot;
block.Derivatives.Data(11)= qdot;
block.Derivatives.Data(12)= rdot;

end 


%% Terminate:
%   Functionality    : Called at the end of simulation for cleanup
%   Required         : Yes
%   C-MEX counterpart: mdlTerminate
%
function Terminate(block)

end 

