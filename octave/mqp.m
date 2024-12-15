% Carregar o arquivo .mat
load('dados.mat');

% Definir as variáveis u e y a partir das colunas de Z
u = Z(:, 2); % Segunda coluna de Z é o sinal de entrada
y = Z(:, 3); % Terceira coluna de Z é o sinal de saída

% Construir a matriz Psi e o vetor y1
Psi = [y(1:end-3), y(2:end-2), y(3:end-1), u(2:end-2), u(3:end-1)];
y1 = y(4:end);

% Definir os pesos para MQP
w1 = 1; 
n1 = 100; % Número de amostras com peso w1
w2 = 0.05; 
n2 = numel(y) - n1 - 3; % Número restante de amostras com peso w2
W = diag([w1 * ones(n1, 1); w2 * ones(n2, 1)]); % Matriz diagonal de pesos

% Estimativa dos parâmetros (MQP)
tetaChapeu = (Psi' * W * Psi) \ (Psi' * W * y1);

% Estimação de yChapeu (dados de estimação)
yChapeu = Psi * tetaChapeu;
csiMQP = y1 - yChapeu; % Resíduos ponderados
JMQP = csiMQP' * W * csiMQP; % Função de custo ponderada

% Estimação de yChapeu1 (infinitos passos à frente)
yChapeu1 = y;
for i = 4:length(yChapeu1)
    yChapeu1(i) = ...
        tetaChapeu(1) * yChapeu1(i-3) + ...
        tetaChapeu(2) * yChapeu1(i-2) + ...
        tetaChapeu(3) * yChapeu1(i-1) + ...
        tetaChapeu(4) * u(i-2) + ...
        tetaChapeu(5) * u(i-1);
end

% Cálculo dos resíduos e função de custo para yChapeu1
csiMQP1 = y1 - yChapeu1(4:end);
JMQP1 = csiMQP1' * W * csiMQP1;

% Plot dos resultados
figure;
plot(0:length(y)-1, y, 'k.-', 'LineWidth', 3, 'MarkerSize', 10); hold on; grid on;
plot(3:length(y)-1, yChapeu, 'r--', 'LineWidth', 3);
plot(0:length(y)-1, yChapeu1, 'g-', 'LineWidth', 3);

% Configuração do gráfico
xlabel('k', 'FontSize', 16);
ylabel('y(k)', 'FontSize', 16);
legend({'Dados de estimação', 'Um passo à frente (MQP)', 'Inf. passos à frente (MQP)'}, ...
       'FontSize', 14, 'Location', 'SouthEast');
title(sprintf('MQP: JMQP = %.2f, JMQP1 = %.2f', JMQP, JMQP1), 'FontSize', 16);
