close all; 
clear;
clc;

% Quadrotor 
W_tot = 4.12;  % N / 420 g

CD = 0.97;
S = 0.01;  % m^2
D = 8*0.0254;  % m propeller diameter 8 inches
A = 4*(pi/4)*D^2;  % m^2 for 4 propellers

% Battery
no_cells = 3;
capacity = 1500;  % mA-hr
E_b = no_cells*3.7*capacity*3600/1000;  % Lecture 5 page 53

% air density 
rho = 1.225;  % kg/m^3
% efficiency 
n_m = 0.75;  % Assuming no motor losses
n_e = 0.85;  % Assuming no ESC losses

%Velocity and Drag
V = 0:1:20;
Drag = 1/2*rho*V.^2*CD*S; %Drag at each velocity p.16/29
alphaD = atan(D./W_tot);

% power required to hover from momentum theory 

for i = 1:length(V)
    Drag(i)= 1/2*rho*V(i)^2*CD*S; %Drag at each velocity p.16/29
    alphaD(i) = atan(Drag(i)/W_tot);
    T = sqrt(W_tot^2+Drag(i)^2);

    A1 = 1;
    A2 = 2 * V(i) * sin(alphaD(i));
    A3 = V(i)^2;
    A4 = 0;
    A5 = -((W_tot^2+Drag(i)^2)/(2*rho*A)^2);
    
   % induced velocity 
   v = roots([A1 A2 A3 A4 A5]);
   for j =(1:4)
        if imag(v(j))== 0 % assume that the imaginary part of x is 0
            if real(v(j))>0 % assume if the real part is positive
                P_tot(i)=T*(v(j)+V(i)*sin(alphaD(i)));
                P_tot_V(i)=P_tot(i)/V(i);
            end
        end
    end
end
% hover endurance from 0th order battery model
% Lecture 7 page 8
[M2,I2]=min(P_tot);
Vmax_endurance = V(I2); %Velocity at max endurance is found at min (P_tot)
P_prop_end = P_tot(I2); %P_tot when V is VmaxEndurance, power per prop
te_max=E_b*n_m*n_e/(P_prop_end) % Max Endurance



% flight range 
[M,I]=min(P_tot_V);
Vmax_range=V(I); %Velocity at max range is found at min (P_tot/V)
P_prop = P_tot(I); %P_tot when V is VmaxRange, power per prop
t_e = E_b*n_m*n_e/(P_prop); 
R_max = t_e*Vmax_range % Max Range





