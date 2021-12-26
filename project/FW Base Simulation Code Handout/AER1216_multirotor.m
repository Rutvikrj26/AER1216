%% AER1216 MULTI_ROTOR CODE
clc;
clear all;

A_roll=[4.2683 -3.1716;4 0];
B_roll=[2;0];
C_roll=[0.7417 0.4405];
D_roll=0;
%
A_pitch=[-3.9784 -2.9796;4 0];
B_pitch=[2;0];
C_pitch=[1.2569 0.6083];
D_pitch=0;
%
A_yaw=-0.0059;
B_yaw=1;
C_yaw=1.2653;
D_yaw=0;
%
A_h=[-5.8200 -3.6046e-6
3.8147e-6 0];
B_h=[1024;0];
C_h=[1.4907e-4 1.3191e3];
D_h=0;
%
A_p2u=-0.665;
B_p2u=2;
C_p2u=-3.0772;
D_p2u=0;
%
A_r2v=-0.4596;
B_r2v=2;
C_r2v=2.3868;
D_r2v=0;
%% 
