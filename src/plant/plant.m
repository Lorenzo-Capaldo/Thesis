clear all
clc
close all
addpath(fileparts(pwd));
z33 = zeros(3,3);

% Sim vars
SIMTIME = 3600;
SIMSTEP = 0.1;
simsize = SIMTIME/SIMSTEP+1;
n_states = 6;
n_obsv = 3;
time = 0:SIMSTEP:SIMTIME;

% x = [Np, Ep, psip, u, v, r, tauN, tauE, tauPSI]
% Ship characteristics
[M, D, Athr, ocean_coeff, wind_coeff] = loadvars();
A = [z33,    eye(3);
     z33, -inv(M)*D];
B = [z33, inv(M)]';
Bv = [z33 z33;
      inv(M), inv(M)];
C = [eye(3,6)];

%wave dynamics
wave_ampl = [0, 0, 0];
wave_length = [70 70 70];

%bias dynamics (current + wind)
% current
V_curr = 2;
beta_c = deg2rad(1); %angle of the current/wind wrt NED
%wind
V_wind = 0;
beta_w = beta_c;
%{
drag_coeff = calc_drag_coeff(angle_of_current, psi);
bias = ocean_curr_coeff .* drag_coeff .* Vc^2; + wind_coeff .* drag_coeff .* Vwind^2;
%}

eig(A)
rank(ctrb(A, B))
rank(obsv(A, C))

x0 = zeros(n_states, 1);
ref = zeros(3,1);
bias = [1e5;0;0];
bias_switch = -1;
sim("plant.slx")

Np = states.Data(:,1); Ep = states.Data(:,2); PSIp = states.Data(:,3);

figure
subplot(311), plot(time, Np)
subplot(312), plot(time, Ep)
subplot(313), plot(time, PSIp*180/pi)


%{
function drag_coeff = calc_drag_coeff(angle_of_wind, heading)
    angle_of_attack = heading-angle_of_wind-180;
    cx = 0.7;
    cy = 0.8; 
    cn = 0.15;
    Cx = -cx*cos(heading);
    Cy = cy*sin(heading);
    Cn = cn*sin(2*heading);
    drag_coeff = [Cx, Cy, Cn]';
end
%}