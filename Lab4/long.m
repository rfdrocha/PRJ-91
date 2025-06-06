clear all
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
phi_0 = -2.08/180*pi;     % rad
theta_0 = 0.59/180*pi;    % rad
psi_0 = -0.02/180*pi;     % rad
m = 2096;          % kg
w0 = 0.32;         % m/sec
u0 = 30.86;        % m/sec
 
A_long = [Xu, Xw, Xq-w0, -g*cos(theta_0);
          Zu, Zw, Zq+u0, -g*cos(phi_0)*sin(theta_0);
          Mu, Mw, Mq, 0;
          0, 0, cos(phi_0), 0];
      
B_long = [X_deltaB; Z_deltaB; M_deltaB; 0];

%% Item A
x0 = [0;0;0;0];

U = @(t)[0.5*2.54*(t>1) - 0.5*2.54*(t>2)];
tspan = 0:0.01:30;
opts = odeset('RelTol',1e-2,'AbsTol',1e-4);

[t,X] = ode45(@(t,X)odefun(t,X,U,A_long,B_long),tspan,x0,opts)

w0 = 0.32;         % m/sec
u0 = 30.86;        % m/sec
theta_0 = 0.59*pi/180;    % rad

u = X(:,1) + u0;
w = X(:,2) + w0;
q = X(:,3);
theta = X(:,4) + theta_0;

figure(1)
plot(t,u*1.94384)
hold on
grid on
xlabel('t (s)');
ylabel('u (kts)');
hold off

figure(2)
plot(t,w*1.94384)
hold on
grid on
xlabel('t (s)');
ylabel('w (kts)');
hold off

figure(3)
plot(t,q*180/pi)
hold on
grid on
xlabel('t (s)');
ylabel('q (deg/s)');
hold off

figure(4)
plot(t,theta*180/pi)
hold on
grid on
xlabel('t (s)');
ylabel('theta (deg)');
hold off


%% Item b
[autovec,autoval] = eig(A_long);
disp("Autovalores:")
autoval = diag(autoval)

%% Item c
omega_n = abs(autoval(2))
omega_d = imag(autoval(2))
zeta = -real(autoval(2))/omega_n
t_2A = -log(2)/(zeta*omega_n)

%% Item d
% Pesquisar requisitos no FAR29

%% Item e
% mostrar conta para obter os coeficientes para um helicoptero qualquer

%% Item f
clear Mu Mw   % liberando as variaveis para que sejam simbolos
syms Mu Mw
A_long_sym = [Xu, Xw, Xq-w0, -g*cos(theta_0);
          Zu, Zw, Zq+u0, -g*cos(phi_0)*sin(theta_0);
          Mu, Mw, Mq, 0;
          0, 0, cos(phi_0), 0];

coefs = charpoly(A_long_sym);
A = coefs(1);
B = coefs(2);
C = coefs(3);
D = coefs(4);
E = coefs(5);

RD4 = B*C*D - A*D*D - B*B*E;
% symvar(RD4); % usado para confirmar a ordem das variaveis na fimplicit 
% Converte para funções numéricas
fA   = matlabFunction(A, 'Vars', [Mu, Mw]);
fB   = matlabFunction(B, 'Vars', [Mu, Mw]);
fC   = matlabFunction(C, 'Vars', [Mu, Mw]);
fD   = matlabFunction(D, 'Vars', [Mu, Mw]);
fE   = matlabFunction(E, 'Vars', [Mu, Mw]);
fRD4 = matlabFunction(RD4, 'Vars', [Mu, Mw]);

% Geração de malha
x_lim_lower = -0.25;
x_lim_upper = 0.5;
y_lim_lower = -0.6;
y_lim_upper = 0.15;
[mu_vals, mw_vals] = meshgrid(linspace(x_lim_lower, x_lim_upper, 400), linspace(y_lim_lower, y_lim_upper, 400));

% Avaliação dos coeficientes
A_vals   = fA(mu_vals, mw_vals);
B_vals   = fB(mu_vals, mw_vals);
C_vals   = fC(mu_vals, mw_vals);
D_vals   = fD(mu_vals, mw_vals);
E_vals   = fE(mu_vals, mw_vals);
RD4_vals = fRD4(mu_vals, mw_vals);

% Máscara: alguma função < 0
mask = (A_vals < 0) | (B_vals < 0) | (C_vals < 0) | ...
       (D_vals < 0) | (E_vals < 0) | (RD4_vals < 0);

% Plot
figure(5); clf;
hold on;

% Região onde algum coeficiente é negativo
contourf(mu_vals, mw_vals, mask, [1 1], 'FaceColor', [0.85 0.85 0.85], 'LineStyle', 'none');

% Fronteiras fimplicit de E e RD4
fimplicit(C, [x_lim_lower x_lim_upper y_lim_lower y_lim_upper], 'b', 'LineWidth', 1.5);
fimplicit(D, [x_lim_lower x_lim_upper y_lim_lower y_lim_upper], 'black', 'LineWidth', 1.5);
fimplicit(E, [x_lim_lower x_lim_upper y_lim_lower y_lim_upper], 'cyan', 'LineWidth', 1.5);
fimplicit(RD4, [x_lim_lower x_lim_upper y_lim_lower y_lim_upper], 'r', 'LineWidth', 1.5);

% Ponto marcado
plot(0.0586, 0.0423, 'ko', 'MarkerFaceColor', 'k');

xlabel('M_u');
ylabel('M_w');
grid on;
legend({'','C = 0','D = 0','E = 0', 'RD4 = 0', 'BO-105C'});
axis equal;
hold off;

%%

% funcao auxiliar para calculo da derivada do vetor X
function dxdt = odefun(t,X,U,A_long,B_long)
    dxdt = A_long*X + B_long*U(t);
end


