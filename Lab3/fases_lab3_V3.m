clear all; clc;

%% DADOS HELICÓPTERO

% Peso máximo de decolagem (kg)
MTOW = 720;

% Combustível (L)
tanque = 140;

% Potência máxima disponível (hp)
P_u_hp = 180;

% Rotor principal R/P
R_RP = 3.85;                   % em m
c_RP = 0.18;                   % em m
pas_RP = 3; 
Omega_RP_rpm = 492;            % em rpm
sigma = pas_RP*c_RP/(pi*R_RP); % solidez do R/P

% Rotor de cauda R/C
R_RC = 0.6;             % em m
c_RC = 0.14;            % em m
pas_RC = 2; 
Omega_RC_rpm = 3100;    % em rpm

% Dimensões gerais do helicóptero
h_helicoptero = 2.427;  % em m

% Área equivalente da fuselagem
A_f_ft = 6; % ft2
A_f_m2 = A_f_ft*0.092903; % m2

% Consumo específico constante (lb/hp.h) - a ideia aqui é que serão
%   consumidos 0,5 libras para cada hora usando 1 hp
Ce = 0.5; % lb/hp.hr

%% VALORES DAS CONSTANTES DAS EQUAÇÕES

% Densidade do combustível (kg/m^3)
rho_combustivel = 800;

% Constantes
C_d0 = 0.012;
ki = 1.15;
n_m = 0.85;

% Escolha da condição de análise: sem vento de proa ou com vento
escolha = 2;    % basicamente mudar entre 1 e 2

if escolha == 1
    vento = 0;
else
    vento = 20; % em kts
end

%% ANÁLISE DO PAIRADO

% Atmosfera
%[~, ~, ~, rho] = atmosisa(0);
rho = 1.14549; % valor atualizado para ISA+20 como solicitado
g = 9.81;

% Helicóptero
m = MTOW;
T = m*g;            % N
A = pi*(R_RP^2);    % m2

% Velocidade induzida do pairado no OGE
V_i0 = sqrt(T/(2*rho*A));  % m/s

% Potência induzida no pairado no OGE
P_i0_OGE = T*V_i0;              % W
P_i0_OGE_hp = P_i0_OGE/745.7;   % hp

% Efeito solo (auxílio gráfico)
h_1_ft = 6;                 % Altura do helicóptero na fase 1 em ft
h_1 = 0.3048*h_1_ft;        % Altura do helicóptero na fase 1 em m
z_RP = h_1 + h_helicoptero;
grafico_abscissa = z_RP/(2*R_RP);
grafico_ordenada = 0.87;    % Valor encontrado traçando no gráfico

Delta_hp_EfeitoSolo = P_i0_OGE_hp*(1 - grafico_ordenada);

% Potência induzida no pairado no IGE
P_i0_IGE_hp = ki*(P_i0_OGE_hp - Delta_hp_EfeitoSolo);


% Potência de perfil das pás
C_p0 = (sigma*C_d0/8)*(1);
Omega_RP = Omega_RP_rpm*(2*pi/60);          % rad/s
P_p0 = C_p0*(rho*A*((Omega_RP*R_RP)^3));    % W
P_p0_hp = P_p0/745.7;                       % hp

% Potência de miscelânea
P_ER_hp = P_p0_hp + P_i0_IGE_hp;
P_misc0_hp = (1 - n_m)*(1/n_m)*P_ER_hp;

% Potência total
P_total_0_hp = P_i0_IGE_hp + P_p0_hp + P_misc0_hp;

% Consumo de combustivel
Consumo_pairado = Ce*P_total_0_hp*(5/60);                               % consumo em libras
Consumo_pairado_litros = Consumo_pairado*0.453592/rho_combustivel*1000; % litros

fprintf('--------------------------------\n');
fprintf('VOO PAIRADO: \n\n');
fprintf('Potencia induzida   - %7.4f hp \n',P_i0_IGE_hp);
fprintf('Potencia de perfil  - %7.4f hp \n',P_p0_hp);
fprintf('Potencia parasita   - %7.4f hp \n',0);
fprintf('Potencia miscelanea - %7.4f hp \n',P_misc0_hp);
fprintf('Potencia de subida  - %7.4f hp \n',0);
fprintf('Potencia de descida - %7.4f hp \n',0);
fprintf('Potencia total      - %7.3f hp \n',P_total_0_hp);
fprintf('Consumo combustivel - %7.3f lt \n',Consumo_pairado_litros);
fprintf('Peso atual - %4.3f kg \n',m);
fprintf('Densidade considerada - %4.3f kg/m^3 \n',rho);
fprintf('--------------------------------\n\n');

%% Subida
m = MTOW - Consumo_pairado*0.453592; % Atualizando peso da aeronave
T = m*g;
rho = 1.0628; % considerado ISA+20 e Zp 2500 ft (ponto médio do segmento)
CT = T/(rho*A*(Omega_RP*R_RP)^2);

v = 1000*0.00508;                   % m/s
v_nivelado = 0;
V_i0_subida = sqrt(T/(2*rho*A));    % m/s

velocidades = linspace(0,140);         % velocidades em kts
velocidades_ms = velocidades*0.514444; % m/s

lambda_i0 = V_i0_subida/(Omega_RP*R_RP);
mus = velocidades_ms/(Omega_RP*R_RP);
lambda_c = v/(Omega_RP*R_RP);
lambda_c_nivelado = v_nivelado/(Omega_RP*R_RP);
lambda_is = zeros(1,100);
lambda_is_nivelado = zeros(1,100);

for i=1:100
    lambda_is(i) = calculate_lambda_i_V2(lambda_i0,mus(i),lambda_c);
    lambda_is_nivelado(i) = calculate_lambda_i_V2(lambda_i0,mus(i),lambda_c_nivelado);
end

% Calculando a Vy
% Calulo de potencia induzida
Coef_Potencias_induzidas_subida = ki*CT^2./(2*sqrt(mus.^2 + (lambda_c + lambda_is).^2));    % Adimensional
Potencias_induzidas_subida = Coef_Potencias_induzidas_subida*rho*A*(Omega_RP*R_RP)^3;       % W
Potencias_induzidas_subida_hp = Potencias_induzidas_subida/745.7;                           % hp

Coef_Potencias_induzidas_subida_nivelado = ki*CT^2./(2*sqrt(mus.^2 + (lambda_c_nivelado + lambda_is_nivelado).^2)); % Adimensional
Potencias_induzidas_subida_nivelado = Coef_Potencias_induzidas_subida_nivelado*rho*A*(Omega_RP*R_RP)^3;             % W
Potencias_induzidas_subida_nivelado_hp = Potencias_induzidas_subida_nivelado/745.7;                                 % hp

% Calculo de potencia de perfil
Coef_Potencias_de_perfil_subida = (sigma*C_d0/8)*(1 + 4.65*mus.^2);                                 % Adimensional
Potencias_de_perfil_subida_hp = Coef_Potencias_de_perfil_subida*(rho*A*(Omega_RP*R_RP)^3)/745.7;    % hp

% Calculo de potencia parasita
Coef_Potencias_parasita_subida = A_f_m2/(2*A)*mus.^3;                                           % Adimensional
Potencias_parasita_subida_hp = Coef_Potencias_parasita_subida*(rho*A*(Omega_RP*R_RP)^3)/745.7;  % hp

% Calculo de potencia de subida
Potencia_subida_hp = T*v/745.7; % hp (já foi definido no roteiro o v)

% Calculo de potencia de miscelanea
P_ER_subida_hp = Potencias_induzidas_subida_hp + Potencias_de_perfil_subida_hp + Potencias_parasita_subida_hp + Potencia_subida_hp;
Potencias_miscelanea_subida_hp = (1 - n_m)*(1/n_m)*P_ER_subida_hp;

P_ER_subida_nivelado_hp = Potencias_induzidas_subida_nivelado_hp + Potencias_de_perfil_subida_hp + Potencias_parasita_subida_hp + Potencia_subida_hp;
Potencia_miscelanea_subida_nivelado_hp = (1 - n_m)*(1/n_m)*P_ER_subida_nivelado_hp;

% Calculo de potencia total
Potencia_total_subida_hp = Potencias_induzidas_subida_hp + Potencias_de_perfil_subida_hp + Potencias_parasita_subida_hp + Potencias_miscelanea_subida_hp + Potencia_subida_hp;

Potencia_total_subida_nivelado_hp = Potencias_induzidas_subida_nivelado_hp + Potencias_de_perfil_subida_hp + Potencias_parasita_subida_hp + Potencia_miscelanea_subida_nivelado_hp + Potencia_subida_hp;

% Verificando velocidade com potencia total mínima:
[Potencia_min_corrigido, idx_min_corrigido] = min(Potencia_total_subida_hp);
[Potencia_min_nivelado, idx_min_nivelado] = min(Potencia_total_subida_nivelado_hp);

% Grafico de potencia de subida
figure(1)
plot(velocidades,Potencia_total_subida_hp,'LineWidth',2);
hold on
plot(velocidades,Potencia_total_subida_nivelado_hp,'--','LineWidth',2);
plot(velocidades,Potencias_miscelanea_subida_hp,'LineWidth',1);
plot(velocidades,Potencias_induzidas_subida_hp,'LineWidth',1);
plot(velocidades,Potencias_de_perfil_subida_hp,'LineWidth',1);
plot(velocidades,Potencias_parasita_subida_hp,'LineWidth',1);
yline(180,'LineWidth',1); % potencia maxima
plot([1 1]*velocidades(idx_min_corrigido), [0 1]*Potencia_min_corrigido, '--b');
plot([1 1]*velocidades(idx_min_nivelado), [0 1]*Potencia_min_nivelado, '--r');
plot(velocidades(idx_min_corrigido), Potencia_min_corrigido, 'ob');
plot(velocidades(idx_min_nivelado), Potencia_min_nivelado, 'or');

grid on
legend({'Pot\^encia Total','Pot\^encia total sem corre\c{c}\~ao','Pot\^encia de miscel\^anea','Pot\^encia induzida','Pot\^encia de perfil','Pot\^encia parasita','Pot\^encia disponivel'},'Location','southwest',Interpreter='latex');
xlabel('Velocidade (kt)',Interpreter='latex');
ylabel('Pot\^encia (hp)',Interpreter='latex');
title('Subida a  1000 ft/min',Interpreter='latex');
set(gca, 'FontSize', 24,'TickLabelInterpreter','latex')
hold off



Consumo_subida = Ce*Potencia_min_corrigido*(5/60); % consumo em libras
Consumo_subida_litros = Consumo_subida*0.453592/rho_combustivel*1000; % litros

fprintf('--------------------------------\n');
fprintf('SUBIDA: \n\n');
fprintf('Potencia induzida   - %7.4f hp \n',Potencias_induzidas_subida_hp(idx_min_corrigido));
fprintf('Potencia de perfil  - %7.4f hp \n',Potencias_de_perfil_subida_hp(idx_min_corrigido));
fprintf('Potencia parasita   - %7.4f hp \n',Potencias_parasita_subida_hp(idx_min_corrigido));
fprintf('Potencia miscelanea - %7.4f hp \n',Potencias_miscelanea_subida_hp(idx_min_corrigido));
fprintf('Potencia de subida  - %7.4f hp \n',Potencia_subida_hp);
fprintf('Potencia de descida - %7.4f hp \n',0);
fprintf('Potencia total      - %7.3f hp \n',Potencia_total_subida_hp(idx_min_corrigido));
fprintf('Consumo combustivel - %7.3f lt \n',Consumo_subida_litros);
fprintf('Vy    - %3.1f kts \n',velocidades(idx_min_corrigido));


% Polar de velocidade
Deltas_P_subida_subida = (P_u_hp - Potencia_total_subida_hp)*745.7;
vs_subida_subida = Deltas_P_subida_subida/(m*g);
vs_subida_subida_kts = vs_subida_subida/0.514444;

Delta_P_descida_subida = -Potencia_total_subida_nivelado_hp*745.7;
vs_descida_subida = Delta_P_descida_subida/(m*g);
vs_descida_subida_kts = vs_descida_subida/0.514444;

% Verificando velocidade com potencia total mínima:
[v_max_subida, idx_max] = max(vs_subida_subida_kts);
[v_min_subida, idx_min] = max(vs_descida_subida_kts);

fprintf('V_vM  - %3.1f kts \n',velocidades(idx_max));
fprintf('V_vm  - %3.1f kts \n',velocidades(idx_min));
fprintf('V_max - %3.1f kts \n',vs_subida_subida_kts(idx_max));
fprintf('V_min - %3.1f kts \n',vs_descida_subida_kts(idx_min));
fprintf('Peso atual - %4.3f kg \n',m);
fprintf('Densidade considerada - %4.3f kg/m^3 \n',rho);
fprintf('--------------------------------\n\n');

% Gráfico da polar
figure(2)
hold on
plot(velocidades,vs_subida_subida_kts,'LineWidth',2)
plot(velocidades,vs_descida_subida_kts,'LineWidth',2)
plot([1 1]*velocidades(idx_max), [0 1]*v_max_subida, '--b');
plot([0 1]*velocidades(idx_max), [1 1]*v_max_subida, '--b');
plot(velocidades(idx_max), v_max_subida, 'ob');
plot(velocidades(idx_max), v_max_subida, 'ob');
plot([1 1]*velocidades(idx_min), [0 1]*v_min_subida, '--r');
plot([0 1]*velocidades(idx_min), [1 1]*v_min_subida, '--r');
plot(velocidades(idx_min), v_min_subida, 'or');
plot([0 1]*velocidades(end), [0 0]*v_min_subida, 'LineWidth', 1,'Color','k')

grid on
xlabel('Velocidade (kt)',Interpreter='latex');
ylabel('v (kt)',Interpreter='latex');
title('Subida com Zp = 2500 ft',Interpreter='latex');
set(gca, 'FontSize', 24,'TickLabelInterpreter','latex')
hold off

%% Nivelado 1
rho = 0.98476; % 5000 ft
m = m - Consumo_subida*0.453592; % Atualizando peso da aeronave
T = m*g;
CT = T/(rho*A*(Omega_RP*R_RP)^2);
V_i0_nivelado_1 = sqrt(T/(2*rho*A)); % m/s

v = 0;

%velocidades = linspace(0,140);         % velocidades em kts
%velocidades_ms = velocidades*0.514444; % m/s

lambda_i0 = V_i0_nivelado_1/(Omega_RP*R_RP);
%mus = velocidades_ms/(Omega_RP*R_RP);
lambda_c = v/(Omega_RP*R_RP);
lambda_is = zeros(1,100);

for i=1:100
    lambda_is(i) = calculate_lambda_i_V2(lambda_i0,mus(i),lambda_c);
end


% Potencia induzida
Coef_Potencias_induzidas_nivelado1 = ki*CT^2./(2*sqrt(mus.^2 + (lambda_c + lambda_is).^2));         % Adimensional
Potencias_induzidas_nivelado_1 = Coef_Potencias_induzidas_nivelado1*rho*A*(Omega_RP*R_RP)^3;        % W
Potencias_induzidas_nivelado_1_hp = Potencias_induzidas_nivelado_1/745.7;                           % hp

% Potencia de perfil
Potencia_perfil_nivelado_1 = rho*A*((Omega_RP*R_RP)^3)*sigma*C_d0/8*(1 + 4.65*mus.^2);  % W
Potencia_perfil_nivelado_1_hp = Potencia_perfil_nivelado_1/745.7;                       % hp

% Potencia parasita
Potencia_parasita_nivelado_1_hp = (A_f_m2/(2*A)*mus.^3)*(rho*A*(Omega_RP*R_RP)^3)/745.7; % hp

% Potencia subida
Potencia_subida_nivelado_1 = v*T;

% Potencia descida
Potencia_descida_nivelado_1 = v*T;

% Potencia miscelanea
P_ER_nivelado_1_hp = Potencias_induzidas_nivelado_1_hp + Potencia_perfil_nivelado_1_hp + Potencia_parasita_nivelado_1_hp; % hp
Potencia_miscelanea_nivelado_1_hp = (1 - n_m)*(1/n_m)*P_ER_nivelado_1_hp; % hp

% Potencia total
Potencia_total_nivelado_1_hp = Potencias_induzidas_nivelado_1_hp + Potencia_perfil_nivelado_1_hp + Potencia_parasita_nivelado_1_hp + Potencia_miscelanea_nivelado_1_hp; % hp

% Calculando a velocidade de maior alcance (depende do vento)


if escolha == 1
    [razao_min, idx_min] = min(Potencia_total_nivelado_1_hp./velocidades);
else
    razao_min = 1e8;        % Apenas um referencial para começar a iteração
    idx_min = 0;
    for i = 1:100
        razao = Potencia_total_nivelado_1_hp(i)/(velocidades(i)-vento);
        if (razao < razao_min) && (razao > 0)
            razao_min = razao;
            idx_min = i;    % Marca o índice correspondente a tangente
        end
    end
end

% Gráfico de potência nivelado 1
figure(3)
plot(velocidades,Potencia_total_nivelado_1_hp,'LineWidth',2);
hold on
plot(velocidades,Potencia_miscelanea_nivelado_1_hp,'LineWidth',1);
plot(velocidades,Potencias_induzidas_nivelado_1_hp,'LineWidth',1);
plot(velocidades,Potencia_perfil_nivelado_1_hp,'LineWidth',1);
plot(velocidades,Potencia_parasita_nivelado_1_hp,'LineWidth',1);
yline(180,'LineWidth',1); % potencia maxima
plot(velocidades(1,1:85)+vento,razao_min*velocidades(1,1:85),'LineWidth',1,'Color','red','LineStyle','--');
plot([1 1]*velocidades(idx_min), [0 1]*Potencia_total_nivelado_1_hp(idx_min), '--r');
plot(velocidades(idx_min), Potencia_total_nivelado_1_hp(idx_min), 'or');

grid on
legend({'Pot\^encia Total','Pot\^encia de miscel\^anea','Pot\^encia induzida','Pot\^encia de perfil','Pot\^encia parasita','Pot\^encia disponivel'},'Location','southwest',Interpreter='latex');
xlabel('Velocidade (kt)',Interpreter='latex');
ylabel('Pot\^encia (hp)',Interpreter='latex');
title('Nivelado na VDM com Zp = 5000 ft',Interpreter='latex');
set(gca, 'FontSize', 24,'TickLabelInterpreter','latex')
hold off

% Consumo nivelado 1
alcance = 300*1852; % m
tempo_nivelado_1 = alcance/(velocidades_ms(idx_min))/60; % tempo em minutos
Consumo_nivelado_1 = Ce*Potencia_total_nivelado_1_hp(idx_min)*(tempo_nivelado_1/60); % consumo em libras
Consumo_nivelado_1_litros = Consumo_nivelado_1*0.453592/rho_combustivel*1000; % litros

fprintf('--------------------------------\n');
fprintf('NIVELADO 1: \n\n');
fprintf('Potencia induzida   - %7.4f hp \n',Potencias_induzidas_nivelado_1_hp(idx_min));
fprintf('Potencia de perfil  - %7.4f hp \n',Potencia_perfil_nivelado_1_hp(idx_min));
fprintf('Potencia parasita   - %7.4f hp \n',Potencia_parasita_nivelado_1_hp(idx_min));
fprintf('Potencia miscelanea - %7.4f hp \n',Potencia_miscelanea_nivelado_1_hp(idx_min));
fprintf('Potencia de subida  - %7.4f hp \n',0);
fprintf('Potencia de descida - %7.4f hp \n',0);
fprintf('Potencia total      - %7.3f hp \n',Potencia_total_nivelado_1_hp(idx_min));
fprintf('Consumo combustivel - %7.3f lt \n',Consumo_nivelado_1_litros);
fprintf('VDM   - %3.1f kts \n',velocidades(idx_min));


% Polar de velocidade
Deltas_P_subida_nivelado_1 = (P_u_hp - Potencia_total_nivelado_1_hp)*745.7;
vs_subida_nivelado_1 = Deltas_P_subida_nivelado_1/(m*g);
vs_subida_nivelado_1_kts = vs_subida_nivelado_1/0.514444;

Delta_P_descida_nivelado_1 = -Potencia_total_nivelado_1_hp*745.7;
vs_descida_nivelado_1 = Delta_P_descida_nivelado_1/(m*g);
vs_descida_nivelado_1_kts = vs_descida_nivelado_1/0.514444;

% Verificando velocidade com potencia total mínima:
[v_max_nivelado_1, idx_max] = max(vs_subida_nivelado_1_kts);
[v_min_nivelado_1, idx_min] = max(vs_descida_nivelado_1_kts);

fprintf('V_vM  - %3.1f kts \n',velocidades(idx_max));
fprintf('V_vm  - %3.1f kts \n',velocidades(idx_min));
fprintf('V_max - %3.1f kts \n',vs_subida_nivelado_1_kts(idx_max));
fprintf('V_min - %3.1f kts \n',vs_descida_nivelado_1_kts(idx_min));
fprintf('Peso atual - %4.3f kg \n',m);
fprintf('Densidade considerada - %4.3f kg/m^3 \n',rho);
fprintf('--------------------------------\n\n');


% Gráfico da polar
figure(4)
hold on
plot(velocidades,vs_subida_nivelado_1_kts,'LineWidth',2)
plot(velocidades,vs_descida_nivelado_1_kts,'LineWidth',2)
plot([1 1]*velocidades(idx_max), [0 1]*v_max_nivelado_1, '--b');
plot([0 1]*velocidades(idx_max), [1 1]*v_max_nivelado_1, '--b');
plot(velocidades(idx_max), v_max_nivelado_1, 'ob');
plot([1 1]*velocidades(idx_min), [0 1]*v_min_nivelado_1, '--r');
plot([0 1]*velocidades(idx_min), [1 1]*v_min_nivelado_1, '--r');
plot(velocidades(idx_min), v_min_nivelado_1, 'or');
plot([0 1]*velocidades(end), [0 0]*v_min_nivelado_1, 'LineWidth', 1,'Color','k')

grid on
xlabel('Velocidade (kt)',Interpreter='latex');
ylabel('v (kt)',Interpreter='latex');
title('Nivelado na VDM com Zp = 5000 ft',Interpreter='latex');
set(gca, 'FontSize', 24,'TickLabelInterpreter','latex')
hold off

%% Nivelado 2
rho = 0.98476; % 5000 ft
m = m - Consumo_nivelado_1*0.453592; % Atualizando peso da aeronave
T = m*g;
CT = T/(rho*A*(Omega_RP*R_RP)^2);
V_i0_nivelado_2 = sqrt(T/(2*rho*A)); % m/s

v = 0;

%velocidades = linspace(0,140);         % velocidades em kts
%velocidades_ms = velocidades*0.514444; % m/s

lambda_i0 = V_i0_nivelado_2/(Omega_RP*R_RP);
%mus = velocidades_ms/(Omega_RP*R_RP);
lambda_c = v/(Omega_RP*R_RP);
lambda_is = zeros(1,100);

for i=1:100
    lambda_is(i) = calculate_lambda_i_V2(lambda_i0,mus(i),lambda_c);
end

% Potencia induzida
%Potencia_induzida_nivelado_2 = ki*(m*g)^2./(2*rho*A*sqrt(velocidades_ms.^2 + V_i0_nivelado_2^2)); % W
%Potencia_induzida_nivelado_2 = Potencia_induzida_nivelado_2/745.7; % hp
Coef_Potencias_induzidas_nivelado2 = ki*CT^2./(2*sqrt(mus.^2 + (lambda_c + lambda_is).^2));     % Adimensional
Potencias_induzidas_nivelado_2 = Coef_Potencias_induzidas_nivelado2*rho*A*(Omega_RP*R_RP)^3;    % W
Potencias_induzidas_nivelado_2_hp = Potencias_induzidas_nivelado_2/745.7;                       % hp

% Potencia de perfil
Potencia_perfil_nivelado_2 = rho*A*((Omega_RP*R_RP)^3)*sigma*C_d0/8*(1 + 4.65*mus.^2);  % W
Potencia_perfil_nivelado_2_hp = Potencia_perfil_nivelado_2/745.7;                       % hp

% Potencia parasita
Potencia_parasita_nivelado_2_hp = (A_f_m2/(2*A)*mus.^3)*(rho*A*(Omega_RP*R_RP)^3)/745.7; % hp

% Potencia subida
Potencia_subida_nivelado_2_hp = v*T;

% Potencia descida
Potencia_descida_nivelado_2_hp = v*T;

% Potencia miscelanea
P_ER_nivelado_2_hp = Potencias_induzidas_nivelado_2_hp + Potencia_perfil_nivelado_2_hp + Potencia_parasita_nivelado_2_hp; % hp
Potencia_miscelanea_nivelado_2_hp = (1 - n_m)*(1/n_m)*P_ER_nivelado_2_hp; % hp

% Potencia total
Potencia_total_nivelado_2_hp = Potencias_induzidas_nivelado_2_hp + Potencia_perfil_nivelado_2_hp + Potencia_parasita_nivelado_2_hp + Potencia_miscelanea_nivelado_2_hp; % hp

% Calculando a velocidade de maior autonomia
[Potencia_min_nivelado_2, idx_min] = min(Potencia_total_nivelado_2_hp);

% Consumo nivelado 2
tempo_nivelado_2 = 20; % tempo em minutos
Consumo_nivelado_2 = Ce*Potencia_total_nivelado_2_hp(idx_min)*(tempo_nivelado_2/60); % consumo em libras
Consumo_nivelado_2_litros = Consumo_nivelado_2*0.453592/rho_combustivel*1000; % litros


% Gráfico de potência nivelado 2
figure(5)
plot(velocidades,Potencia_total_nivelado_2_hp,'LineWidth',2);
hold on
plot(velocidades,Potencia_miscelanea_nivelado_2_hp,'LineWidth',1);
plot(velocidades,Potencias_induzidas_nivelado_2_hp,'LineWidth',1);
plot(velocidades,Potencia_perfil_nivelado_2_hp,'LineWidth',1);
plot(velocidades,Potencia_parasita_nivelado_2_hp,'LineWidth',1);
yline(180,'LineWidth',1); % potencia maxima

plot([1 1]*velocidades(idx_min), [0 1]*Potencia_min_nivelado_2, '--r');
plot(velocidades(idx_min), Potencia_min_nivelado_2, 'or');

grid on
legend({'Pot\^encia Total','Pot\^encia de miscel\^anea','Pot\^encia induzida','Pot\^encia de perfil','Pot\^encia parasita','Pot\^encia disponivel'},'Location','southwest',Interpreter='latex');
xlabel('Velocidade (kt)',Interpreter='latex');
ylabel('Pot\^encia (hp)',Interpreter='latex');
title('Nivelado na VAM com Zp = 5000 ft',Interpreter='latex');
set(gca, 'FontSize', 24,'TickLabelInterpreter','latex')
hold off

% Imprimindo resultados nivelado 2

fprintf('--------------------------------\n');
fprintf('NIVELADO 2: \n\n');
fprintf('Potencia induzida   - %7.4f hp \n',Potencias_induzidas_nivelado_2_hp(idx_min));
fprintf('Potencia de perfil  - %7.4f hp \n',Potencia_perfil_nivelado_2_hp(idx_min));
fprintf('Potencia parasita   - %7.4f hp \n',Potencia_parasita_nivelado_2_hp(idx_min));
fprintf('Potencia miscelanea - %7.4f hp \n',Potencia_miscelanea_nivelado_2_hp(idx_min));
fprintf('Potencia de subida  - %7.4f hp \n',0);
fprintf('Potencia de descida - %7.4f hp \n',0);
fprintf('Potencia total      - %7.3f hp \n',Potencia_total_nivelado_2_hp(idx_min));
fprintf('Consumo combustivel - %7.3f lt \n',Consumo_nivelado_2_litros);
fprintf('VAM   - %3.1f kts \n',velocidades(idx_min));


% Polar de velocidade
Deltas_P_subida_nivelado_2 = (P_u_hp - Potencia_total_nivelado_2_hp)*745.7;
vs_subida_nivelado_2 = Deltas_P_subida_nivelado_2/(m*g);
vs_subida_nivelado_2_kts = vs_subida_nivelado_2/0.514444;

Delta_P_descida_nivelado_2 = -Potencia_total_nivelado_2_hp*745.7;
vs_descida_nivelado_2 = Delta_P_descida_nivelado_2/(m*g);
vs_descida_nivelado_2_kts = vs_descida_nivelado_2/0.514444;

% Verificando velocidade com potencia total mínima:
[v_max_nivelado_2, idx_max] = max(vs_subida_nivelado_2_kts);
[v_min_nivelado_2, idx_min] = max(vs_descida_nivelado_2_kts);

fprintf('V_vM  - %3.1f kts \n',velocidades(idx_max));
fprintf('V_vm  - %3.1f kts \n',velocidades(idx_min));
fprintf('V_max - %3.1f kts \n',vs_subida_nivelado_2_kts(idx_max));
fprintf('V_min - %3.1f kts \n',vs_descida_nivelado_2_kts(idx_min));
fprintf('Peso atual - %4.3f kg \n',m);
fprintf('Densidade considerada - %4.3f kg/m^3 \n',rho);
fprintf('--------------------------------\n\n');

% Gráfico da polar
figure(6)
hold on
plot(velocidades,vs_subida_nivelado_2_kts,'LineWidth',2)
plot(velocidades,vs_descida_nivelado_2_kts,'LineWidth',2)
plot([1 1]*velocidades(idx_max), [0 1]*v_max_nivelado_2, '--b');
plot([0 1]*velocidades(idx_max), [1 1]*v_max_nivelado_2, '--b');
plot(velocidades(idx_max), v_max_nivelado_2, 'ob');
plot([1 1]*velocidades(idx_min), [0 1]*v_min_nivelado_2, '--r');
plot([0 1]*velocidades(idx_min), [1 1]*v_min_nivelado_2, '--r');
plot(velocidades(idx_min), v_min_nivelado_2, 'or');
plot([0 1]*velocidades(end), [0 0]*v_min_nivelado_2, 'LineWidth', 1,'Color','k')

grid on
xlabel('Velocidade (kt)',Interpreter='latex');
ylabel('v (kt)',Interpreter='latex');
title('Nivelado na VAM com Zp = 5000 ft',Interpreter='latex');
set(gca, 'FontSize', 24,'TickLabelInterpreter','latex')
hold off

%% Descida na Vy

rho = 1.0628; % considerado ISA+20 e Zp 2500 ft (ponto médio do segmento)
m = m - Consumo_nivelado_2*0.453592; % Atualizando peso da aeronave
T = m*g;
CT = T/(rho*A*(Omega_RP*R_RP)^2);


v = -1000*0.00508; % m/s (negativo pois é descida)
V_i0_descida = sqrt(T/(2*rho*A)); % m/s


velocidades = linspace(0,140);         % velocidades em kts
velocidades_ms = velocidades*0.514444; % m/s

lambda_i0 = V_i0_descida/(Omega_RP*R_RP);
mus = velocidades_ms/(Omega_RP*R_RP);
lambda_c = 0;       % pois não é subida, logo, não há compensação de AUMENTO de m_ponto
%lambda_d = v/(Omega_RP*R_RP); nem precisa desse termo porque calculo
%P_descida direto por Pot = T*V
lambda_is = zeros(1,100);

for i=1:100
    lambda_is(i) = calculate_lambda_i_V2(lambda_i0,mus(i),lambda_c);
end

% Calculando a Vy
% Calulo de potencia induzida
Coef_Potencias_induzidas_descida = ki*CT^2./(2*sqrt(mus.^2 + (lambda_c + lambda_is).^2));   % Adimensional
Potencias_induzidas_descida = Coef_Potencias_induzidas_descida*rho*A*(Omega_RP*R_RP)^3;     % W
Potencias_induzidas_descida_hp = Potencias_induzidas_descida/745.7;                         % hp


% Calculo de potencia de perfil
Coef_Potencia_de_perfil_descida = (sigma*C_d0/8)*(1 + 4.65*mus.^2);                              % Adimensional
Potencia_de_perfil_descida_hp = Coef_Potencia_de_perfil_descida*(rho*A*(Omega_RP*R_RP)^3)/745.7; % hp

% Calculo de potencia parasita
Coef_Potencia_parasita_descida = A_f_m2/(2*A)*mus.^3;                                           % Adimensional
Potencia_parasita_descida_hp = Coef_Potencia_parasita_descida*(rho*A*(Omega_RP*R_RP)^3)/745.7;  % hp

% Calculo de potencia de subida
Potencia_subida = 0;

% Calculo de potencia de descida
Potencia_descida_hp = T*v/745.7;    % hp

% Calculo de potencia de miscelanea
P_ER_descida_hp = Potencias_induzidas_descida_hp + Potencia_de_perfil_descida_hp + Potencia_parasita_descida_hp + Potencia_descida_hp;
Potencia_miscelanea_descida_hp = (1 - n_m)*(1/n_m)*P_ER_descida_hp;

% Calculo de potencia total
Potencia_total_descida_hp = Potencias_induzidas_descida_hp + Potencia_de_perfil_descida_hp + Potencia_parasita_descida_hp + Potencia_miscelanea_descida_hp + Potencia_descida_hp;

% Verificando velocidade com potencia total mínima:
[Potencia_min, idx_min] = min(Potencia_total_descida_hp);

Consumo_descida = Ce*Potencia_min*(5/60); % consumo em libras
Consumo_descida_litros = Consumo_descida*0.453592/rho_combustivel*1000; % litros

% Grafico de potencia de subida
figure(7)
plot(velocidades,Potencia_total_descida_hp,'LineWidth',2);
hold on
plot(velocidades,Potencia_miscelanea_descida_hp,'LineWidth',1);
plot(velocidades,Potencias_induzidas_descida_hp,'LineWidth',1);
plot(velocidades,Potencia_de_perfil_descida_hp,'LineWidth',1);
plot(velocidades,Potencia_parasita_descida_hp,'LineWidth',1);
plot(velocidades,Potencia_descida_hp*ones(1, 100),'LineWidth',1);
yline(180,'LineWidth',1); % potencia maxima
plot([1 1]*velocidades(idx_min), [0 1]*Potencia_min, '--r');
plot(velocidades(idx_min), Potencia_min, 'or');

grid on
legend({'Pot\^encia Total','Pot\^encia de miscel\^anea','Pot\^encia induzida','Pot\^encia de perfil','Pot\^encia parasita','Pot\^encia de descida','Pot\^encia disponivel'},'Location','southwest',Interpreter='latex');
xlabel('Velocidade (kt)',Interpreter='latex');
ylabel('Pot\^encia (hp)',Interpreter='latex');
title('Descida a 1000 ft/min',Interpreter='latex');
set(gca, 'FontSize', 20,'TickLabelInterpreter','latex')
hold off

fprintf('--------------------------------\n');
fprintf('DESCIDA: \n\n');
fprintf('Potencia induzida   - %7.4f hp \n',Potencias_induzidas_descida_hp(idx_min));
fprintf('Potencia de perfil  - %7.4f hp \n',Potencia_de_perfil_descida_hp(idx_min));
fprintf('Potencia parasita   - %7.4f hp \n',Potencia_parasita_descida_hp(idx_min));
fprintf('Potencia miscelanea - %7.4f hp \n',Potencia_miscelanea_descida_hp(idx_min));
fprintf('Potencia de subida  - %7.4f hp \n',0);
fprintf('Potencia de descida - %7.4f hp \n',Potencia_descida_hp);
fprintf('Potencia total      - %7.3f hp \n',Potencia_total_descida_hp(idx_min));
fprintf('Consumo combustivel - %7.3f lt \n',Consumo_descida_litros);
fprintf('Vy    - %3.1f kts \n',velocidades(idx_min));


% Polar de velocidade
Deltas_P_subida_descida = (P_u_hp - Potencia_total_descida_hp)*745.7;
vs_subida_descida = Deltas_P_subida_descida/(m*g);
vs_subida_descida_kts = vs_subida_descida/0.514444;

Delta_P_descida_descida = -Potencia_total_descida_hp*745.7;
vs_descida_descida = Delta_P_descida_descida/(m*g);
vs_descida_descida_kts = vs_descida_descida/0.514444;

% Verificando velocidade com potencia total mínima:
[v_max_descida, idx_max] = max(vs_subida_descida_kts);
[v_min_descida, idx_min] = max(vs_descida_descida_kts);

fprintf('V_vM  - %3.1f kts \n',velocidades(idx_max));
fprintf('V_vm  - %3.1f kts \n',velocidades(idx_min));
fprintf('V_max - %3.1f kts \n',vs_subida_descida_kts(idx_max));
fprintf('V_min - %3.1f kts \n',vs_descida_descida_kts(idx_min));
fprintf('Peso atual - %4.3f kg \n',m);
fprintf('Densidade considerada - %4.3f kg/m^3 \n',rho);
fprintf('--------------------------------\n\n');

% Gráfico da polar
figure(8)
hold on
plot(velocidades,vs_subida_descida_kts,'LineWidth',2)
plot(velocidades,vs_descida_descida_kts,'LineWidth',2)
plot([1 1]*velocidades(idx_max), [0 1]*v_max_descida, '--b');
plot([0 1]*velocidades(idx_max), [1 1]*v_max_descida, '--b');
plot(velocidades(idx_max), v_max_descida, 'ob');
plot([1 1]*velocidades(idx_min), [0 1]*v_min_descida, '--r');
plot([0 1]*velocidades(idx_min), [1 1]*v_min_descida, '--r');
plot(velocidades(idx_min), v_min_descida, 'or');
plot([0 1]*velocidades(end), [0 0]*v_min_descida, 'LineWidth', 1,'Color','k')

grid on
xlabel('Velocidade (kt)',Interpreter='latex');
ylabel('v (kt)',Interpreter='latex');
title('Descida com Zp = 2500 ft',Interpreter='latex');
set(gca, 'FontSize', 24,'TickLabelInterpreter','latex')
hold off


%% Pairado IGE final
rho = 1.14549;
m = m - Consumo_descida*0.453592; % Atualizando peso da aeronave
T = m*g;

% Velocidade induzida do pairado no OGE
V_i0 = sqrt(T/(2*rho*A));  % m/s

% Potência induzida no pairado no OGE
P_i0_OGE = T*V_i0;   % W
P_i0_OGE_hp = P_i0_OGE/745.7; % hp

% Correcao Ground Effect
grafico_ordenada = 0.87;    % Valor encontrado traçando no gráfico
Delta_hp_EfeitoSolo = P_i0_OGE_hp*(1 - grafico_ordenada);

% Potência induzida no pairado no IGE
P_i0_IGE_hp = ki*(P_i0_OGE_hp - Delta_hp_EfeitoSolo);

% Potência de perfil das pás
C_p0 = (sigma*C_d0/8)*(1);
Omega_RP = Omega_RP_rpm*(2*pi/60);
P_p0 = C_p0*(rho*A*((Omega_RP*R_RP)^3)); 
P_p0_hp = P_p0/745.7; 

% Potência de miscelânea
P_ER = P_p0_hp + P_i0_IGE_hp;
P_misc0 = (1 - n_m)*(1/n_m)*P_ER;

% Potência total
P_total_0 = P_i0_IGE_hp + P_p0_hp + P_misc0;

% Consumo de combustivel
Consumo_pairado = Ce*P_total_0*(5/60); % consumo em libras
Consumo_pairado_litros = Consumo_pairado*0.453592/rho_combustivel*1000; % litros

fprintf('--------------------------------\n');
fprintf('VOO PAIRADO 2: \n\n');
fprintf('Potencia induzida   - %7.4f hp \n',P_i0_IGE_hp);
fprintf('Potencia de perfil  - %7.4f hp \n',P_p0_hp);
fprintf('Potencia parasita   - %7.4f hp \n',0);
fprintf('Potencia miscelanea - %7.4f hp \n',P_misc0);
fprintf('Potencia de subida  - %7.4f hp \n',0);
fprintf('Potencia de descida - %7.4f hp \n',0);
fprintf('Potencia total      - %7.3f hp \n',P_total_0);
fprintf('Consumo combustivel - %7.3f lt \n',Consumo_pairado_litros);
fprintf('Peso atual - %4.3f kg \n',m);
fprintf('Densidade considerada - %4.3f kg/m^3 \n',rho);
fprintf('--------------------------------\n\n');