clear all
close all
clc

echo on
% SIMcruise User editable script for simulation of the cruise ship 
% % under feedback control. Based on T. I. Fossen's MSS scripts.
%
% Calls:       cruise.m, Lcruise.m and euler2.m
%
% Author:      L. Capaldo
% Date:        2023-09-24
% Revisions: 

echo off
disp('Simulating cruise.m under PID control with psi_ref = 5 deg ...')

% SIMPARAMS
t_f = 600;   % final simulation time (sec)
h   = 0.1;   % sample time (sec)
N = round(t_f/h);                       % number of samples

% Memory allocation
% x_vec = [Np Ep psip u v r bn bE bPSI]
% y_vec = [N E psi]
% out   = [time Np Ep psip u v r bn bE bPSI tau_X tau_Y tau_N N E psi udot vdot rdot U]
n_states  = 9;
n_inputs  = 3;
n_outputs = 3;
n_acce    = 3;
x    = zeros(n_states, N);
y    = zeros(n_outputs, N);
out  = zeros(1+n_states+n_inputs+n_outputs+n_acce+1, N);          
%tau_vec = zeros(3, N+1);

% Initial values
x(:, 1) = [0 0 0*pi/180 0 0 0 0 0 0]'; 
x(:, 2) = [0 0 0*pi/180 0 0 0 0 0 0]'; 
tau =[0, 0, 0]';
%tau_vec(1, 1:end) = ones(size(1, 200))*400000;
out(:,1) = [1; x(:,1); tau; y(:,1); [0;0;0]; sqrt(x(4, 1)^2 + x(5, 1)^2)];
out(:,2) = [2; x(:,2); tau; y(:,2); [0;0;0]; sqrt(x(4, 2)^2 + x(5, 2)^2)];
tau_prev = tau;

% MAIN LOOP
for i=3:N

    time = i*h;                   % simulation time in seconds

    % control system
    % tau = -Kp * ( ssa(psi-psi_ref) + Td * r);   % PD controller
    pos_ref = [100, 100, 0]';                         % desired position
    tau = control(x(1:6, i-1), y(:,i-1), y(:,i-2), pos_ref, 100, 1, 100);

    % ship model
    [xdot, U] = cruise(x(:, i-1), tau);          % ship model
    
    % Euler integration
    x(:, i) = x(:, i-1) + xdot * h;

    % Convert to NED for the output
    psip = x(3, i); % = psi
    rotmat = [cos(psip), -sin(psip), 0; 
              sin(psip),  cos(psip), 0;
                      0,          0, 1];
    y(:, i) = inv(rotmat') * x(1:3, i); % convert to normal NED

    % store data for presentation
    out(:, i) = [time; x(:, i); tau; y(:, i); xdot(4:6); U]';

end

plot_all(out, pos_ref)

function tau = control(x_1, y_1, y_2, pos_ref, Kp, Kd, Ki)
    Kp = eye(3) * Kp;      % controller P gain
    Kd = eye(3) * Kd;
    Ki = eye(3) * Ki;

    delta_eta = y_1 - pos_ref;
    psip = y_1(3); % = psi
    rotmat = [cos(psip), -sin(psip), 0; 
              sin(psip),  cos(psip), 0;
                      0,          0, 1];

    tau = - rotmat' * Kp * delta_eta - rotmat' * Kd * rotmat * x_1(4:6) - rotmat' * Ki * (y_2-pos_ref);
end

function plot_all(xout, pos_ref)
    % time-series
    % x_vec = [Np Ep psip u v r bn bE bPSI]
    % y_vec = [N E psi]
    % out   = [time Np Ep psip u v r bn bE bPSI tau_X tau_Y tau_N N E psi udot vdot rdot U]
    t         = xout(1,:);
    xp        = xout(2,:); %x
    yp        = xout(3,:); %y
    psip      = xout(4,:) * 180/pi;
    u         = xout(5,:); 
    v         = xout(6,:);          
    r         = xout(7,:) * 180/pi;  
    b_N       = xout(8,:);
    b_E       = xout(9,:);
    b_psi     = xout(10,:);
    tau_X     = xout(11,:);
    tau_Y     = xout(12,:);
    tau_N     = xout(13,:);
    x         = xout(14,:);
    y         = xout(15,:);
    psi       = xout(16,:);
    udot      = xout(17,:);
    vdot      = xout(18,:);
    rdot      = xout(19,:);
    U         = xout(20,:);

    % Plots
    %{
    figure
    plot(y,x,'r')
    hold on
    grid,axis('equal'),xlabel('East'),ylabel('North'),title('Ship position')
    plot(y(1), x(1),'bx');
    plot(y(end), x(end),'ro');
    hold off
    set(findall(gcf,'type','line'),'linewidth',2)
    set(findall(gcf,'type','text'),'FontSize',14)
    set(findall(gcf,'type','legend'),'FontSize',14)
    %}
    
    figure
    subplot(311),plot(t,x,'r'),xlabel('time (s)'),title('position x (m)'),grid
    hold on, plot(t, pos_ref(1)*ones(size(t)), 'b--'), hold off
    subplot(312),plot(t,y,'r'),xlabel('time (s)'),title('position y (m)'),grid
    hold on, plot(t, pos_ref(2)*ones(size(t)), 'b--'), hold off
    subplot(313),plot(t,psi,'r'),xlabel('time (s)'),title('yaw angle \psi (deg)'),grid
    hold on, plot(t, pos_ref(3)*ones(size(t)), 'b--'), hold off
    set(findall(gcf,'type','line'),'linewidth',2)
    set(findall(gcf,'type','text'),'FontSize',14)
    set(findall(gcf,'type','legend'),'FontSize',14)
    
    %{
    figure
    subplot(221), plot(t, u, 'r'), xlabel('time (s)'), title('u (m/s)'),grid
    subplot(222), plot(t, v, 'r'), xlabel('time (s)'), title('v (m/s)'),grid
    subplot(223), plot(t, r, 'r'), xlabel('time (s)'), title('r (deg/s)'),grid
    subplot(224), plot(t, U, 'r'), xlabel('time (s)'), title('Speed (m/s)'),grid
    set(findall(gcf,'type','line'),'linewidth',2)
    set(findall(gcf,'type','text'),'FontSize',14)
    set(findall(gcf,'type','legend'),'FontSize',14)
    
    figure
    subplot(311), plot(t, udot, 'r'), xlabel('time (s)'), title('$\dot u$ ', 'Interpreter','latex'),grid
    subplot(312), plot(t, vdot, 'r'), xlabel('time (s)'), title('$\dot v$', 'Interpreter','latex'),grid
    subplot(313), plot(t, rdot, 'r'), xlabel('time (s)'), title('$\dot r$', 'Interpreter','latex'),grid
    set(findall(gcf,'type','line'),'linewidth',2)
    set(findall(gcf,'type','text'),'FontSize',14)
    set(findall(gcf,'type','legend'),'FontSize',14)
    
    figure
    subplot(311), plot(t, b_N, 'r'), xlabel('time (s)'), title('bN'),grid
    subplot(312), plot(t, b_E, 'r'), xlabel('time (s)'), title('bE'),grid
    subplot(313), plot(t, b_psi, 'r'), xlabel('time (s)'), title('bPSI'),grid
    set(findall(gcf,'type','line'),'linewidth',2)
    set(findall(gcf,'type','text'),'FontSize',14)
    set(findall(gcf,'type','legend'),'FontSize',14)
    %}
    
    
    figure
    subplot(311), plot(t, tau_X, 'r'), xlabel('time (s)'), title('x thrust (N)'),grid
    subplot(312), plot(t, tau_Y, 'r'), xlabel('time (s)'), title('y thrust (N)'),grid
    subplot(313), plot(t, tau_N, 'r'), xlabel('time (s)'), title('\psi thrust (N)'),grid
    set(findall(gcf,'type','line'),'linewidth',2)
    set(findall(gcf,'type','text'),'FontSize',14)
    set(findall(gcf,'type','legend'),'FontSize',14)
    
end
