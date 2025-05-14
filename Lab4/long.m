% Constantes
Xu = -0.0338;      % 1/sec
Xw = 0.0311;       % 1/sec
Xq = 0.6384;       % m/sec-rad
X_deltaB = 0.0844; % m/sec2-cm
Zu = -0.0564;      % 1/sec
Zw = -0.7886;      % 1/sec
Zq = 0.0565;       % m/sec-rad
Z_deltaB = 0.2157; % m/sec2-cm
Mu = 0.0586;       % rad/m-sec
Mw = 0.0423;       % rad/m-sec
Mq = -3.6151;      % 1/sec
M_deltaB = -0.3922;% rad/sec2-cm
g = 9.81;          % m/s2
phi_0 = -2.08;     % deg
theta_0 = 0.59;    % deg
psi_0 = -0.02;     % deg
m = 2096;          % kg
w0 = 0.32;         % m/sec
u0 = 30.86;        % m/sec
 
A_long = [Xu, Xw, Xq-w0, -g*cos(theta_0);
          Zu, Zw, Zq+u0, -g*cos(phi_0)*sin(theta_0);
          Mu, Mw, Mq, 0;
          0, 0, cos(phi_0), 0];
      
B_long = [X_deltaB; Z_deltaB; M_deltaB; 0];


tspan = [0 30];
x0 = [0;0;0;0];

