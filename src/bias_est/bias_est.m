clear all
clc
close all
addpath(fileparts(pwd)+"\loadvars.m");
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

Q = diag([1e3 1e3 1e3 10 10 10]);
R = diag([.1 .1 .1]);
K = lqr(A, B, Q, R);
cl_eig = eig(A-B*K);
%
% Bias Estimation
n_aug = 3;
Aa = [A, Bv(:,1:3); zeros(3, 6), z33];
Ba = [B;z33]; Ca = [C, z33];
rank(obsv(Aa, Ca))

lambda_obs = [5.*cl_eig' -.8 -.2 -0.3];
L = place(Aa', Ca', lambda_obs)';

xhat_0 = zeros(n_states+n_aug,1);
x0 = zeros(n_states,1);
u = zeros(n_inputs,1);
bias = [1e4;1e3;0];
sim("biasest.slx")

states_hat = logsout{1}.Values.Data';
bias_sim = squeeze(logsout{2}.Values.Data);
states = logsout{3}.Values.Data';

figure, sgtitle("Pos vs posest")
subplot(311), plot(time, states(1,:)), hold on, plot(time, states_hat(1,:), "r--")
subplot(312), plot(time, states(2,:)), hold on, plot(time, states_hat(2,:), "r--")
subplot(313), plot(time, states(3,:)), hold on, plot(time, states_hat(3,:), "r--")

figure, sgtitle("Bias Estimation")
subplot(311), plot(time, bias_sim(1,:)), hold on, plot(time, states_hat(7,:), "r--")
subplot(312), plot(time, bias_sim(2,:)), hold on, plot(time, states_hat(8,:), "r--")
subplot(313), plot(time, bias_sim(3,:)), hold on, plot(time, states_hat(9,:), "r--")

%%
% Observer + Controller
% First we need to find the controllable subspace
n_aug = 3;
Aa = [A, Bv(:,1:3); zeros(3, 9)];
Ba = [B;zeros(3,3)]; Ca = [C, zeros(3,3)];

dim_Anc = size(Aa, 1) - rank(ctrb(Aa, Ba));
disp('System uncontrollable subspace has dimension '), disp(dim_Anc) 
disp('Controllable subspace decomposition')
[Aac,Bac,Cac,Q,~] = ctrbf(Aa,Ba,Ca);
Ac = Aac(dim_Anc+1:end,dim_Anc+1:end);
Bc = Bac(dim_Anc+1:end,:);
Cc = Cac(:, dim_Anc+1:end);

lambda_ctrl_des = [-1 -1.2 -0.6 -0.8 -0.5 -0.4]';
Kt1_1 = place(Ac, Bc, lambda_ctrl_des);
Kt2_1 = [0 0 0;0 0 0;0 0 0]; 
Kt2_2 = [0.5 0.5 0.5 ;0.5 0.5 0.5; 0.5 0.5 0.5]; 
Kt1 = [Kt2_1 Kt1_1];
Kt2 = [Kt2_2 Kt1_1];
K1 = Kt1*Q;
K2 = Kt2*Q;
eig(Ac-Bc*Kt1_1)
xhat_0 = zeros(n_states+n_aug,1);
x0 = zeros(n_states,1);
ref = zeros(n_outputs, 1);
bias = [1e4;1e3;0];
sim("biasest_with_control.slx")

states_hat = logsout{1}.Values.Data';
bias_sim = squeeze(logsout{2}.Values.Data);
states = logsout{3}.Values.Data';

figure, sgtitle("Pos vs posest")
subplot(311), plot(time, states(1,:)), hold on, plot(time, states_hat(1,:), "r--")
subplot(312), plot(time, states(2,:)), hold on, plot(time, states_hat(2,:), "r--")
subplot(313), plot(time, states(3,:)), hold on, plot(time, states_hat(3,:), "r--")

figure, sgtitle("Bias Estimation")
subplot(311), plot(time, bias_sim(1,:)), hold on, plot(time, states_hat(7,:), "r--")
subplot(312), plot(time, bias_sim(2,:)), hold on, plot(time, states_hat(8,:), "r--")
subplot(313), plot(time, bias_sim(3,:)), hold on, plot(time, states_hat(9,:), "r--")