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
n_inputs = 3;
n_outputs = 3;
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
bias = [1e4;1e3;0];
bias_switch = +1; % -1 for constant bias, +1 for dynamic bias
ref = zeros(3,1);
sim("plant.slx")

Np = logsout{1}.Values.Data(:,1); Ep = logsout{1}.Values.Data(:,2); PSIp = logsout{1}.Values.Data(:,3);
N = logsout{1}.Values.Data(:,1); E = logsout{1}.Values.Data(:,2); PSI = logsout{1}.Values.Data(:,3);

figure
subplot(311), plot(time, Np), hold on, plot(time, N, "r--")
subplot(312), plot(time, Ep), hold on, plot(time, E, "r--")
subplot(313), plot(time, PSIp), hold on, plot(time, PSI, "r--")

%%
n_augm = 3;
Ai = [A zeros(6,3);
      -C zeros(3,3)];
Bi = [B; zeros(3,3)];
Bvi = [Bv; zeros(3, 6)];
Ci = [C zeros(3,3)];
rank(ctrb(Ai, Bi))

Q = diag([1e5, 1e5, 1e5, 0.01, 0.01, 0.01, 10, 10, 10]);
R = diag([0.01, 0.01, 0.01]);
K = lqr(Ai, Bi, Q, R);
Klqr = K(:,1:n_states);
Ki = -K(:,end-n_augm+1:end);

ref = [2 1 deg2rad(1)];
bias = [1e4;1e3;0];
bias_switch = -1;
x0 = zeros(n_states, 1);
x0(1:3) = [1, 1, 0];
x0int = zeros(n_augm ,1);
sim("lqr_int_sim.slx")

states_NED = squeeze(logsout{1}.Values.Data);

figure, sgtitle("Controlled Position")
subplot(311), plot(time, states_NED(1,:))
subplot(312), plot(time, states_NED(2,:))
subplot(313), plot(time, states_NED(3,:))