%% SAMPLE FILE FOR AER1216 ---- Q4 b)
%
clc;
clear all;
%
%%  Basic Inputs
load para.mat
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
k = 1 / (pi*e*AR);

C_T=0.78; % We need to redefine


%  (2) flight conditions reference point ****

Va_trim=15;
Va = 28.35;
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
Xdelta_t = 0;

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
%%
%Lateral 
% X_latr=[v p r phi psi]';
% U_latr=[delta_a delta_r]';

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
%    Blatr(5,:)=[0 0];
% X_latr_dot=Alatr*X_latr + Blatr*U_latr;
Ta=20.2;
Tr=20.2;
u_v=0;
u_psi=180;
A_latr_new=[Alatr Blatr zeros(5,2);zeros(1,5) -Ta 0 0 0;zeros(1,5) 0 -Tr 0 0;-1 zeros(1,8);zeros(1,4) -1 zeros(1,4)];
%
B_latr_new=[zeros(5,2);Ta 0;0 Tr;zeros(2,2)];
E_latr_new=[zeros(7,2);1 0;0 1];
E_latr_new=E_latr_new*[u_v u_psi]';
%%
C_latr_new=eye(9);
D_latr_new=zeros(9,2);
F_latr_new=zeros(9,2);
F_latr_new(8:9,:)=eye(2);
F_latr_new=F_latr_new*[u_v u_psi]';
H_latr_new=[0 0 0 0 1 0 0 0 0];
%%
state_7=ss(A_latr_new,B_latr_new,C_latr_new,D_latr_new);
Qq=0.0001*eye(9)+H_latr_new'*H_latr_new; %working
%%
R1=100;%working
[Y1,pp,e2]=lqr(state_7,Qq,R1);
% [X1,S]=linsolve(C_latr_new',Y1');
K_latr_new=Y1;
%%
Kx_latr=K_latr_new(:,1:5);
Kdel_latr=K_latr_new(:,6:7);
Keps_latr=K_latr_new(:,8:9);
%closed loop to verify results
A3=A_latr_new-B_latr_new*K_latr_new*C_latr_new;
B3=E_latr_new-B_latr_new*K_latr_new*F_latr_new;

% The new closed loop system is %%MIGHT NOT NEED
sys_new2=ss(A3,B3,C_latr_new,F_latr_new);
figure(3)
t=0:0.01:1000;
[mm1,l1,x3]=step(sys_new2,t);
del_a=(180/pi)*[0 0 0 0 0 1 0 0 0]*x3';
psi=[0 0 0 0 1 0 0 0 0]*x3';
v=[1 0 0 0 0 0 0 0 0]*x3';
phi=[0 0 0 1 0 0 0 0 0]*x3';
%% plotting the step command
plot(t,v,t,del_a,t,psi,t,phi) %
grid
legend('v','delta a','psi','phi');
xlabel('Time (sec)');
title('Step Response');
