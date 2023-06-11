close all; clear; clc;

%% Given Parameter
P.m = 1200;
P.L = 2.5;
P.Caf = 5500;
P.Car = 4500;
P.Iz = 1800;
P.R = 210;
P.wdis = 1.5;
P.g = 9.81;
P.L_s = 500;

%% Derived Parameter
P.lf = P.L*(1 - P.wdis/(1 + P.wdis));
P.lr = P.L*P.wdis/(1 + P.wdis);


%% Set Parameter
P.tou = 0.500;
P.Kp = 1.1667;
P.Ki = 0.4815;

%% State Space
state0 = zeros(1,14);
state0(1) = 1;
state0(3) = 10/180*pi;
state0(12) = 10;
tspan = [0, 60]; 
opts = odeset('MaxStep',0.01);

dstatedt = @(t, state) statespace(t, state, P);
[t, state] = ode45 (dstatedt, tspan, state0, opts);

%% plot
% e1, e2
figure(1);
subplot(211);
plot(t, state(:,1), '-', 'linewidth', 1.5); hold on;
grid on; grid minor;
ylabel('Lateral Displacement [m]', 'interpreter', 'latex', 'fontsize', 12);
title('e1', 'interpreter', 'latex', 'fontsize', 12);
subplot(212);
plot(t, state(:,3),'-','linewidth',1.5); hold on;
grid on; grid minor;
ylabel('Angle [rad]', 'interpreter', 'latex', 'fontsize', 12);
xlabel('Time [sec]', 'interpreter', 'latex', 'fontsize', 12);
title('e2', 'interpreter', 'latex', 'fontsize', 12);

%Trajectory
figure(2);
plot(state(:,6), state(:,7), 'LineWidth', 2); hold on;
plot(state(:,8), state(:,9), '--r','LineWidth', 2); hold on;
grid on; grid minor;
legend('Actual trajectory', 'Desired Trajectory', 'interpreter', 'latex', 'fontsize', 12);
xlabel('x [m]', 'interpreter', 'latex', 'fontsize', 12);
ylabel('y [m]', 'interpreter', 'latex', 'fontsize', 12);
title('Inertial Trajectory');
xlim([ 0 1000]);
ylim([-100 1000]);

% Cruise Control
for i = 1:length(t)
   v_des(i) = x_dot_des(t(i));
end
figure(3)
plot(t, state(:, 12), '-b', 'LineWidth', 1.5); hold on;
plot(t, v_des, '--r', 'LineWidth', 1.5); hold on;
grid on; grid minor;
legend('$v$', '$v_{des}$', 'interpreter', 'latex', 'fontsize', 12);
xlabel('Time [sec]', 'interpreter', 'latex', 'fontsize', 12);
ylabel('Velocity [m/s]', 'interpreter', 'latex', 'fontsize', 12);
title('Cruise Control', 'fontsize', 12);
ylim([0 50]);

%% ODE function of state
function dstatedt = statespace(t, state, P)
  dstatedt = zeros(14,1);
  %state = [e1, e1', e2, e2', ]
  Vx = state(12);

  %% state_ e1, e1', e2, e2', psi_des, X, Y, Xc, Yc, e1_i, e2_i, Vx, Vx', Vx''
  A = [ 0                                       1                                 0                                           0;
        0             -(2*P.Caf+2*P.Car)/(P.m*Vx)             (2*P.Caf+2*P.Car)/P.m       (-2*P.Caf*P.lf+2*P.Car*P.lr)/(P.m*Vx);
        0                                       0                                 0                                           1;
        0  -(2*P.Caf*P.lf-2*P.Car*P.lr)/(P.Iz*Vx)  (2*P.Caf*P.lf-2*P.Car*P.lr)/P.Iz  -(2*P.Caf*P.lf^2+2*P.Car*P.lr^2)/(P.Iz*Vx)];        

  B = [              0                                             0    0;
             2*P.Caf/P.m  -(2*P.Caf*P.lf-2*P.Car*P.lr)/(P.m*Vx) - Vx  P.g;
                       0                                           0    0;
       2*P.Caf*P.lf/P.Iz  -(2*P.Caf*P.lf^2+2*P.Car*P.lr^2)/(P.Iz*Vx)    0];

  %% delta
  delta = str_ang(state);

  %% phi_d_des
  if 0 < state(8) && state(8) < P.L_s
    phi_d_des = 0;
  elseif P.R <= state(8) && state(9) < P.R
    phi_d_des = Vx/P.R;
  else 
    phi_d_des = 0;
  end

  %% input & calculate
  x = state(1:4);

  u = [    delta;
       phi_d_des;
               0]; 
  
  dstatedt(1:4) = A*x + B*u;

  %% state_ psi_des, X, Y
  dstatedt(5) = phi_d_des;
  dstatedt(6) = Vx*cos(state(5)) - (state(2)*sin(state(3)+state(5)) + state(1)*cos(state(3)+state(5))*(state(4) + phi_d_des));
  dstatedt(7) = Vx*sin(state(5)) + (state(2)*cos(state(3)+state(5)) - state(1)*sin(state(3)+state(5))*(state(4) + phi_d_des));
  
  %% state_ Xc, Yc, e1_i, e2_i
  dstatedt(8) = Vx*cos(state(5));
  dstatedt(9) = Vx*sin(state(5));
  dstatedt(10) = state(1);
  dstatedt(11) = state(3);

  %% Cruise control state_ Vx, Vx', Vx''
  Av = [          0           1        0;
                  0           0        1;
        -P.Ki/P.tou -P.Kp/P.tou -1/P.tou];
  Bv = [         0; 
                 0; 
        P.Ki/P.tou];
  dstatedt(12:14) = Av * state(12:14) + Bv * x_dot_des(t);

end