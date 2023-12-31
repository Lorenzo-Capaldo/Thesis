function [xdot, U] = cruise_fossen(x,tau)
% [xdot, U] = supply(x,tau) returns the speed and the time derivative 
%  xdot = A*x + B*tau of the state vector: x = [ x y psi u v r ]'  for a 
% supply vessel length L = 76 m.
%
% The model is only valid for staionkeeping (dynamic positioning) and 
% low speed U = sqrt( u^2 + v^2 ).
%
% u:   surge velocity                    (m/s)     
% v:   sway velocity                     (m/s)
% r:   yaw velocity                      (rad/s)
% x    position in the x-direction       (m)
% y:   position in the y-direction       (m)
% psi: yaw angle                         (rad)
%
% tau   = [X, Y, N]' control force/moment
%
% Reference: Fossen, T. I., S. I. Sagatun and A. J. Sorensen (1996)
%            Identification of Dynamically Positioned Ships
%            Journal of Control Engineering Practice CEP-4(3):369-376
%
% Author:     Thor I. Fossen
% Date:       12 July 2002
% Revisions:  24 February 2004 - Included missing mass in the Bis transformation
%             12 October 2011 - Corrected T and Tinv, which were switched 
%             27 May 2019 - Added U as ouput
%             31 May 2019 - Included the rotation matrix in yaw
%             22 Mar 2023 - Corrected wrong assignmnet of states

% Normalization variables
L    =  292.5;           % length of ship (m)
g    =  9.8;            % acceleration of gravity (m/s^2)
mass = 85619e3;          % mass (kg)

T    = diag([1 1 1/L]);
Tinv = diag([1 1 L]);

% Model matricses
Mbis = [1.1274         0          0
             0    1.8902    -0.0744
             0   -0.0744     0.1278];

Dbis = [0.0358        0        0
             0        0.1183  -0.0124
             0       -0.0041   0.0308];
 
 psi = x(3);
 R = [ cos(psi) -sin(psi) 0
       sin(psi)  cos(psi) 0
             0         0  1 ];
  
 M = mass * Tinv^2 * (T * Mbis * Tinv);
 D = mass * Tinv^2 * (sqrt(g/L) * T * Dbis * Tinv);

 A = [ zeros(3,3)         R
       zeros(3,3) -inv(M)*D ];
 
 B = [ zeros(3,3)
       inv(M) ];
 A
 B
 % Dimensional state derivative
 xdot = A * x + B * tau;
 
 % speed
  U = sqrt( x(4)^2 + x(5)^2 );