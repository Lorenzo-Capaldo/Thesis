clear all
clc
close all

% Sim vars
SIMTIME = 3600;
SIMSTEP = 0.1;
simsize = SIMTIME/SIMSTEP+1;
n_contr = 9;
n_obsv = 9;
time = 0:SIMSTEP:SIMTIME;

z33 = zeros(3,3);
% x = [Np, Ep, psip, u, v, r, tauN, tauE, tauPSI]
% Ship characteristics
[M, D, Athr] = loadvars();
A = [z33,    eye(3),    z33;
     z33, -inv(M)*D, inv(M);
     z33,       z33,   Athr];
B = [z33, z33, -Athr]';
Bv = [z33, inv(M), z33]';
Bvw = Bv;
wave_ampl = [2, 2, 0];
wave_length = [70 70 70];
C = [eye(9,9)];
eig(A);
rank(ctrb(A, B));
rank(obsv(A, C));

x0 = zeros(n_contr, 1);
bias = [0, 0, 0];
ref = zeros(3,1);
sim("plant.slx")

Np = y.Data(:,1); Ep = y.Data(:,2); PSIp = y.Data(:,3);

figure
subplot(311), plot(time, Np)
subplot(312), plot(time, Ep)
subplot(313), plot(time, PSIp)

%%

A = [z33,    eye(3),    z33;
     z33, -inv(M)*D, inv(M);
     z33,       z33,   Athr];
B = [z33, z33, -Athr]';
Bv = [z33, inv(M), z33]';
Bvw = Bv;
wave_ampl = [2, 2, 0].*0;
wave_length = [70 70 70];
C = eye(3,9);
n_augm = 3;

Ai = [A zeros(9,3);
      -C zeros(3,3)];
Bi = [B; zeros(3,3)];
N = [zeros(9,3); eye(3)];
Bvi = [Bv; zeros(3, 3)];
Ci = [C zeros(3,3)];
rank(ctrb(Ai, Bi))

Q = diag([1e5, 1e5, 1e5, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 10, 10, 10]);
R = diag([0.01, 0.01, 0.01]);
K = lqr(Ai, Bi, Q, R);
Klqr = K(:,1:n_contr);
Ki = -K(:,end-n_augm+1:end);

ref = [100 0 0];
bias = [5, 1, 0];
x0 = zeros(n_contr, 1);
x0(1:3) = [1, 1, 0];
x0int = zeros(n_augm ,1);
sim("lqr_int.slx")

Np = y.Data(:,1); Ep = y.Data(:,2); PSIp = y.Data(:,3);
tauN = states.Data(:,7); tauE = states.Data(:,8); tauPSI = states.Data(:,9); 

figure
subplot(311), plot(time, Np)
subplot(312), plot(time, Ep)
subplot(313), plot(time, PSIp)
figure, sgtitle("taus")
subplot(311), plot(time, tauN)
subplot(312), plot(time, tauE)
subplot(313), plot(time, tauPSI)