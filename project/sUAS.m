%% 1216 Project Q2

clear all
clc

%  (1) aircraft configuration parameters ****
g = 9.81; 			% gravity (m/s^2)
m = 13.5; 	    	% mass (kg)
W = m*g;   % weight (N)
Ixx = 0.8244;   % moment of inertia (kg m^2)
Iyy = 1.135;   % 
Izz = 1.759;   %
Ixz = 0.1204;
Izx = Ixz;
J = [ Ixx   0 -Ixz
        0 Iyy    0
     -Izx   0  Izz];

S = 0.55; 			% wing plan area (m^2)
c_bar = 0.18994; 	    % mean chord (m) 
b = 2.8956; 			% wing span (m) 
AR = b / c_bar;
epsilon = 0.1592;     % oswald factor, CHECK
k = 1 / (pi*epsilon*AR);
C_D_0 = 0.03;
C_m_0 = -0.02338;
% ht = 2.84; % height (m)

%  (2) flight conditions reference point ****
%% DOUBT
Ue = 67; 		% steady flight velocity (m/s): 220 fps 
%% DOUBT ^^
% alt = 1500;     % altitude (m) @ stdatm 5000ft
%                 % or run [T,P,rho] = stdatm(alt);
% [T,P,rho] = stdatm(alt);
rho = 1.0582; 	% density (kg/m^3) 
%% CHECK ^^
% speed_of_sound = sqrt(1.4*287*T);
speed_of_sound=343; %m/s
Mach = Ue / speed_of_sound;

%  (3) steady states
theta_e = 0;                        % pitch angle (rad)
Cw = W / (0.5 * rho * S * Ue^2);    % non-dimensional term of weight \hat{mg}
CLe = Cw * cos(theta_e); 	        % initial C_L lift coefficients 
% Cze = -0.9530;                    % given by original data, not used here
CDe = 0.032; 	                    % initial C_D lift coefficients (given)
Cze = - CLe;
Cxe = Cw * sin(theta_e);
CTe = CDe + Cxe;
Cme = 0.0;


%   (5) derivatives 

% % (4.1) the following inertia values are reverse-engineered ...
% iyy = Iyy / ( rho * S * (c_bar / 2)^3);     % non-dimensional term of moment of inertia i_{yy}
% ixx = Ixx / ( rho * S * (b / 2)^3);         % non-dimensional term of moment of inertia i_{xx}
% izz = Izz / ( rho * S * (b / 2)^3);         % non-dimensional term of moment of inertia i_{xx}
% ixz = Ixz / ( rho * S * (b / 2)^3);         % non-dimensional term of moment of inertia i_{xx}
% 
% % (4.2) mass
% muc = m / (0.5 * rho * S * c_bar);          % non-dimensional term of mass \mu for longitudinal calculations
% mub = m / (0.5 * rho * S * b);              % non-dimensional term of mass \mu for lateral calculations
% 
Cxu  = 0.03;          % CTu = -0.096; CDu = 0; Cxu = CTu - CDu;
Cxa  = 0.30; 		% C_{x_\alpha}
Cxq  =  0.0;
Cxdu =  0.0;			% C_{x_{\dot u}}
Cxda =  0.0; 			% C_{x_{\dot \alpha}}
Cxdq =  0.0;

Czu  = 0.28; %Cz =-CL
Cza  = 3.45;
Czq  = 0;
Czdu =  0.0;
Czda =  0;
Czdq =  0.0;

Cmu  =-0.02338  ;
Cma  = -0.38;
Cmq  = -3.6 ;
Cmdu =  0.0;
Cmda = -7.27;
Cmdq =  0.0;

Cxdelta_e = 0; %C_{x_{\delta_e}}
CLdelta_e = -0.36;   
Czdelta_e = -CLdelta_e;
Cmdelta_e = -0.5;
Cxdelta_p = 0;  % this and the following delta_p are arbitrary
Czdelta_p = 0;
Cmdelta_p = 0;

Cyb = -0.98;          % C_{y_\beta}
Cyp = 0;
Cyr = 0;
Cydb = 0;               % C_{y_{\dot \beta}}
Cydp = 0;
Cydr = 0;

Clb = -0.12;
Clp = -0.26;
Clr =  0.14;
Cldb = 0;               % C_{l_{\dot \beta}}
Cldp = 0;
Cldr = 0;

Cnb = 0.25;
Cnp = 0.022;
Cnr = -0.35;
Cndb = 0;               % C_{n_{\dot \beta}}
Cndp = 0;
Cndr = 0;

Cydelta_a = 0;          %C_{y_{\delta_a}}
Cydelta_r = -0.17;
Cldelta_a = 0.08;
Cldelta_r = 0.105;
Cndelta_a = 0.06;
Cndelta_r = -0.032;

% =========  dimensional derivatives ================
q = 0.5*rho*Ue^2;
qs = rho * Ue * S / 2;
qs2 = rho * Ue * Ue * S / 2;

% ( Longitudinal terms
Xu  = (Cxu+2*Cxe) * qs;
Xw  = Cxa * qs;
Xq  = Cxq * qs * c_bar / 2;
Xdw = Cxda * rho * S * c_bar / 4;  
Xdu = Cxdu * rho * S * c_bar / 4;
Xdq = Cxdq * rho * S * c_bar^2 / 8;

Zu  = (Czu + 2*Cze) * qs;
Zw  = Cza * qs;
Zq  = Czq * qs * c_bar / 2;
Zdu =  Czdu * rho * S * c_bar / 4;
Zdw = Czda * rho * S * c_bar / 4;
Zdq =  Czdq * rho * S * c_bar^2 / 8;

Mu  = (Cmu + 2*Cme) * qs * c_bar;
Mw  = Cma * qs * c_bar;
Mq  = Cmq * rho * Ue * S * c_bar * c_bar / 4;
Mdu = Cmdu * rho * S * c_bar^2 / 4;
Mdw = Cmda * rho * S * c_bar^2 / 4;
Mdq = Cmdq * rho * S * c_bar^3 / 8;

Xdelta_e = Cxdelta_e * qs * Ue;
Zdelta_e = Czdelta_e * qs * Ue;
Mdelta_e = Cmdelta_e * qs * Ue * c_bar;

Xdelta_p = Cxdelta_p * qs * Ue;
Zdelta_p = Czdelta_p * qs * Ue;
Mdelta_p = Cmdelta_p * qs * Ue * c_bar;

%  Lateral terms
Yv = Cyb * qs;
Yp = Cyp * qs * b / 2;
Yr = Cyr * qs * b / 2;
Ydv = Cydb * rho * S * b / 4;
Ydp = Cydp * rho * S * b^2 / 8;
Ydr = Cydr * rho * S * b^2 / 8;

Lv = Clb * qs * b;
Lp = Clp * qs * b^2 / 2;
Lr = Clr * qs * b^2 / 2;
Ldv = Cldb * rho * S * b^2 / 4;
Ldp = Cldp * rho * S * b^3 / 8;
Ldr = Cldr * rho * S * b^3 / 8;

Nv = Cnb * qs * b;
Np = Cnp * qs * b^2 / 2;
Nr = Cnr * qs * b^2 / 2;
Ndv = Cndb * rho * S * b^2 / 4;
Ndp = Cndp * rho * S * b^3 / 8;
Ndr = Cndr * rho * S * b^3 / 8;

Ydelta_a = Cydelta_a * qs2;
Ydelta_r = Cydelta_r * qs2;
Ldelta_a = Cldelta_a * qs2 * b;
Ldelta_r = Cldelta_r * qs2 * b;
Ndelta_a = Cndelta_a * qs2 * b;
Ndelta_r = Cndelta_r * qs2 * b;

%  Moment of Inertia data
Ixx0 = (Ixx*Izz - Ixz^2) / Izz;
Izz0 = (Ixx*Izz - Ixz^2) / Ixx;
Ixz0 = Ixz / (Ixx*Izz - Ixz^2);

%% Q2
syms u v w p q r theta phi del_e del_p del_r del_a
% Longitudinal
Xw_dot=Xdw  ;
Xu_dot=Xdu ;
Xq_dot=Xdq ;

Zu_dot= Zdu ;
Zw_dot=Zdw ;
Zq_dot=Zdq ;

Mu_dot= Mdu ;
Mw_dot= Mdw ;
Mq_dot= Mdq ;

X_del_e= Xdelta_e ;
Z_del_e= Zdelta_e ;
M_del_e=Mdelta_e ;

X_del_p=Xdelta_p ;
Z_del_p=Zdelta_p ;
M_del_p=Mdelta_p ;

%  Lateral terms
Yv_dot =Ydv ;
Yp_dot= Ydp ;
Yr_dot= Ydr ;

Lv_dot = Ldv ;
Lp_dot= Ldp;
Lr_dot= Ldr ;

Nv_dot =  Ndv ;
Np_dot = Ndp ;
Nr_dot = Ndr ;

Y_del_a=Ydelta_a;
Y_del_r= Ydelta_r ;
L_del_a= Ldelta_a ;
L_del_r=Ldelta_r ;
N_del_a= Ndelta_a ;
N_del_r = Ndelta_r ;

%
Iy=Iyy;
Ix=Ixx;
Iz=Izz;
%%
%Longitudinal 
X_long=[u w q theta]';
U_long=[del_e del_p]';

Along=[Xu/m Xw/m Xq/m -g*cos(theta_e);Zu/(m-Zw_dot) Zw/(m-Zw_dot) (Zq+m*Ue)/(m-Zw_dot) -W*sin(theta_e)/(m-Zw_dot);(1/Iy)*(Mu+(Zu*Mw_dot/(m-Zw_dot))) (1/Iy)*(Mw+(Zw*Mw_dot/(m-Zw_dot))) (1/Iy)*(Mq+((Zq+m*Ue)*Mw_dot/(m-Zw_dot))) (-1/Iy)*(W*sin(theta_e)*Mw_dot/(m-Zw_dot));0 0 1 0];

Blong=[X_del_e/m X_del_p/m;Z_del_e/(m-Zw_dot) Z_del_p/(m-Zw_dot);(1/Iy)*(M_del_e+(Z_del_e*Mw_dot/(m-Zw_dot))) (1/Iy)*(M_del_p+(Z_del_p*Mw_dot/(m-Zw_dot)));0 0];



%Lateral 
X_latr=[v p r phi]';
U_latr=[del_a del_r]';

Ixx_d=Iz/(Ix*Iz-Ixz^2);
Izz_d=Ix/(Ix*Iz-Ixz^2);
Ixz_d=Ixz/(Ix*Iz-Ixz^2);

Alatr=[Yv/m Yp/m (Yr/m)-Ue g*cos(theta_e);Ixx_d*Lv+Ixz_d*Nv Ixx_d*Lp+Ixz_d*Np Ixx_d*Lr+Ixz_d*Nr 0;Izz_d*Nv+Ixz_d*Lv Izz_d*Np+Ixz_d*Lp Izz_d*Nr+Ixz_d*Lr 0;0 1 tan(theta_e) 0];

Blatr=[Y_del_a/m Y_del_r/m; Ixx_d*L_del_a+Ixz_d*N_del_a Ixx_d*L_del_r+Ixz_d*N_del_r;Izz_d*N_del_a+Ixz_d*L_del_a Izz_d*N_del_r+Ixz_d*L_del_r;0 0];


%%

%% Q3 Speed control through delta_e
%Using PID control design 
%% Q3
s=tf('s');
x_1= (inv(s*eye(4) - Along))*Blong; %
%Selected control channel is speed control through delta_t
G=tf(x_1.Numerator(1,1),x_1.Denominator(1,1)); %The transfer function for speed control through delta_t is
figure(1)
% margin(G);
step(G)
title('Open Loop Step Response')
ylabel('Speed (m/s)'); 
grid on
stepinfo(G)

% Design specifications were :
% Classical design approach used was PID control
Kp=0.1; %
Ki=0.1; %
Kd=0.1; %
s=tf('s');
C=Kp+Ki/s+Kd*s;
%
figure(2)
sys=feedback(C*G,1);
step(sys)
title(['Closed Loop Step Response', 'Kp= ',num2str(Kp),' Ki=',num2str(Ki),'Kd=',num2str(Kd)])
ylabel('Speed (m/s)'); 
grid on
stepinfo(sys)

