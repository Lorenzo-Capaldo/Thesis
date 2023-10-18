function [Mtot, D, Athrust] = loadvars()
% nu = [u, v, r]';
% eta = [N, E, psi];
% x = [Np, Ep, psip, u, v, r, bN, bE, bPSI]';

%Physical characteristics
m = 85619e3;
Iz = 47.934e6;
rho = 1e3;
L = 292.5;
xg = L/3;
b = L;
a = 32.2;

M = [m,    0,    0;
     0,    m, m*xg;
     0, m*xg,   Iz];

% Added mass characteristics
Xud = 1.4e6;
Yvd = 6.8e6;
Yrd = 3.9e7;
Nrd = 4.2e9;

Madd = [-Xud,    0,    0;
           0, -Yvd, -Yrd;
           0, -Yrd, -Nrd];
Mtot = M + Madd;

% Time constants
Tsurge = 150;
Tsway  = 200;
Tyaw   = 180;

% Damping characteristics
A11_0 = 1.4e6;
A22_0 = 7.2e6;
A66_0 = 4.2e9;

B11v   = (m+A11_0)/Tsurge;
B22v   = (m+A22_0)/Tsway;
B66v   = (m+A66_0)/Tyaw;

Xu = -B11v;
Yv = -B22v;
Nr = -B66v;
Yr = 2.9e7;
Nv = 9.6e6;

%Xu = -1.1e6;
%Yv = -3.6e6;
%Yr = 2.9e7;
%Nv = 9.6e6;
%Nr = -5.4912e9;

D = -[Xu,   0,   0;
       0,  Yv,  Yr;
       0,  Nv,  Nr];

Athrust = -[1/Tsurge,       0,      0;
                   0, 1/Tsway,      0;
                   0,       0, 1/Tyaw];
end