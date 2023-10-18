function [Mtot, D, Athrust, ocean_currents_coeffs] = loadvars()
% nu = [u, v, r]';
% eta = [N, E, psi];
% x = [Np, Ep, psip, u, v, r, bN, bE, bPSI]';

%Physical characteristics
m = 1480e3;
rho = 1e3;
L = 90;
xg = L/3;
b = L;
a = 13.28;
draft = 4;
Iz = 1/12 * m * (a^2+b^2); %moment of inertia of a rectangle

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
Tsurge = 70;
Tsway  = 100;
Tyaw   = 150;

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

ocean_currents_coeffs = 1/2*rho.*[b*draft; a*draft; a*draft*b];
end