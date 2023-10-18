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
ref = zeros(n_outputs,1);
bias = [1e5;0;0];
bias_switch = -1;
sim("plant.slx")

%%
% il controllore da solo non riesce ad eliminare il disturbo ma almeno
% stabilizza il sistema
Q = diag([1e5, 1e5, 1e5, 0.1, 0.1, 0.1]);
R = diag([1, 1, 1]);
K = lqr(A, B, Q, R);
Ak = A-B*K;
cl_poles = eig(Ak);
x = zeros(n_states, simsize);
u = zeros(n_inputs, simsize);
x(1,1) = 10;
ref = zeros(n_states,1)';
for i = 2:simsize
    u(:,i-1) = -K*(x(:,i-1)-ref');
    xdot = A*x(:,i-1) + B * u(:,i-1) + Bv(:,1:3) * bias;
    x(:,i) = x(:,i-1) + SIMSTEP*xdot;
end
u(:,i) = -K*(x(:,i)-ref');

figure, sgtitle("Pos Con")
subplot(311), plot(time, x(1,:))
subplot(312), plot(time, x(2,:))
subplot(313), plot(time, x(3,:))
figure, sgtitle("Inputs Con")
subplot(311), plot(time, u(1,:))
subplot(312), plot(time, u(2,:))
subplot(313), plot(time, u(3,:))
%%

% Observer
n_aug = 3;
Aa = [A, Bv(:,1:3); zeros(3, 6), z33];
Ba = [B;z33]; Ca = [C, z33];

lambda_obs = [7.*cl_poles' -1 -2 -3];
L = place(Aa', Ca', lambda_obs)';
obs_cl_eig = eig(Aa-L*Ca)
xhat_0 = zeros(n_states+n_aug,1);
x0 = zeros(n_states,1);
u = zeros(n_inputs,1);
ref = zeros(n_outputs,1);
sim("biasest_distobsv.slx")

states_hat = squeeze(logsout{1}.Values.Data);
bias_sim = squeeze(logsout{2}.Values.Data);

figure, sgtitle("Bias vs Bias Estimation")
subplot(311), plot(time, bias_sim(1,:)), hold on, plot(time, states_hat(7,:), "r--")
subplot(312), plot(time, bias_sim(2,:)), hold on, plot(time, states_hat(8,:), "r--")
subplot(313), plot(time, bias_sim(3,:)), hold on, plot(time, states_hat(9,:), "r--")

%%
% Observer + Controller
% First we need to find the controllable subspace

dim_Anc = size(Aa, 1) - rank(ctrb(Aa, Ba));
disp('System uncontrollable subspace has dimension '), disp(dim_Anc) 
disp('Controllable subspace decomposition')
[Aac,Bac,Cac,Q,~] = ctrbf(Aa,Ba,Ca);
Ac = Aac(dim_Anc+1:end,dim_Anc+1:end);
Bc = Bac(dim_Anc+1:end,:);
Cc = Cac(:, dim_Anc+1:end);

Qlqr = diag([1e5, 1e5, 1e5, 1, 1, 1]);
Rlqr = diag([0.1, 0.1, 0.1]);
Kt1_1 = lqr(Ac, Bc, Qlqr, Rlqr);
Kt2_1 = [0 0 0;0 0 0;0 0 0]; 
Kt2_2 = [0.5 0.5 0.5;0.5 0.5 0.5;0.5 0.5 0.5]; 
Kt1 = [Kt2_1 Kt1_1];
Kt2 = [Kt2_2 Kt1_1];
K1 = Kt1*Q;
K2 = Kt2*Q;

xhat_0 = zeros(n_states+n_aug,1);
x0 = zeros(n_states,1);
ref = zeros(n_outputs,1);
ref = [2;0;0];
sim("lqr_dist_obsv_sim.slx")

states_hat = squeeze(logsout{1}.Values.Data);
bias_sim = squeeze(logsout{2}.Values.Data);
states = logsout{3}.Values.Data';

figure, sgtitle("Pos vs posest")
subplot(311), plot(time, states(1,:)), hold on, plot(time, states_hat(1,:), "r--")
subplot(312), plot(time, states(2,:)), hold on, plot(time, states_hat(2,:), "r--")
subplot(313), plot(time, states(3,:)), hold on, plot(time, states_hat(3,:), "r--")

figure, sgtitle("Bias vs Bias Estimation")
subplot(311), plot(time, bias_sim(1,:)), hold on, plot(time, states_hat(7,:), "r--")
subplot(312), plot(time, bias_sim(2,:)), hold on, plot(time, states_hat(8,:), "r--")
subplot(313), plot(time, bias_sim(3,:)), hold on, plot(time, states_hat(9,:), "r--")