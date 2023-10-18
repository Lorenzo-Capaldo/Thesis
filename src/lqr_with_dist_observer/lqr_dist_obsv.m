clear all
clc
close all

% Sim vars
SIMTIME = 3600;
SIMSTEP = 0.1;
simsize = SIMTIME/SIMSTEP+1;
time = 0:SIMSTEP:SIMTIME;

z33 = zeros(3,3);
n_states = 6;
n_inputs = 3;
n_outputs = 3;
% x = [Np, Ep, psip, u, v, r, tauN, tauE, tauPSI]
% Ship characteristics
[M, D, Athr, ocean_currents_coeffs] = loadvars_smol_boat();
A = [z33,    eye(3);
     z33, -inv(M)*D];
B = [z33, inv(M)]';
Bv = [z33, inv(M)]';
C = [eye(n_outputs,n_states)];

eig(A)
rank(ctrb(A, B))
rank(obsv(A, C))

Vc = 2;
angle_attack = 180;
Cx_at_aa = [0.5, 0, 0];
bias_X = ocean_currents_coeffs(1) * Cx_at_aa(1) * Vc^2;
bias = [bias_X;0;0];

x0 = zeros(n_states, 1);
ref = zeros(n_inputs,1);
u = [0,0,0]';
sim("plant.slx")

Np = y.Data(:,1); Ep = y.Data(:,2); PSIp = y.Data(:,3);
velN = states.Data(:,4); velE = states.Data(:,5); valPSI = states.Data(:,6);

figure
subplot(311), plot(time, Np)
subplot(312), plot(time, Ep)
subplot(313), plot(time, PSIp)


% Observer
n_aug = 3;
Aa = [A, Bv; zeros(3, 6), z33];
Ba = [B;z33]; Ca = [C, z33];

lambda_obs = [3.*cl_poles' -1 -2 -3];
L = place(Aa', Ca', lambda_obs)';

xhat_0 = zeros(n_states+n_aug,1);
x0 = zeros(n_states,1);
u = zeros(n_inputs,1);
sim("biasest.slx")

states = logsout{2}.Values.Data';
states_hat = logsout{1}.Values.Data';

figure, sgtitle("Pos vs posest")
subplot(311), plot(time, states(1,:)), hold on, plot(time, states_hat(1,:), "r--")
subplot(312), plot(time, states(2,:)), hold on, plot(time, states_hat(2,:), "r--")
subplot(313), plot(time, states(3,:)), hold on, plot(time, states_hat(3,:), "r--")

figure, sgtitle("Bias Estimation")
subplot(311), plot(time, states_hat(7,:))
subplot(312), plot(time, states_hat(8,:))
subplot(313), plot(time, states_hat(9,:))

%%
% Observer + Controller
% First we need to find the controllable subspace
n_aug = 3;
Aa = [A, Bv; zeros(3, 6), z33];
Ba = [B;z33]; Ca = [C, z33];

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
Kt2_2 = [0.5 0.5 0.5;0.5 0.5 0.5;0.5 0.5 0.5]; 
Kt1 = [Kt2_1 Kt1_1];
Kt2 = [Kt2_2 Kt1_1];
K1 = Kt1*Q;
K2 = Kt2*Q;
eig(Ac-Bc*Kt1_1)
xhat_0 = zeros(n_states+n_aug,1);
x0 = zeros(n_states,1);
ref = zeros(n_outputs, 1);
sim("fsfc_bias_est.slx")

states = logsout{2}.Values.Data';
states_hat = logsout{1}.Values.Data';

figure, sgtitle("Pos vs posest")
subplot(311), plot(time, states(1,:)), hold on, plot(time, states_hat(1,:), "r--")
subplot(312), plot(time, states(2,:)), hold on, plot(time, states_hat(2,:), "r--")
subplot(313), plot(time, states(3,:)), hold on, plot(time, states_hat(3,:), "r--")

figure, sgtitle("Bias Estimation")
subplot(311), plot(time, states_hat(7,:))
subplot(312), plot(time, states_hat(8,:))
subplot(313), plot(time, states_hat(9,:))