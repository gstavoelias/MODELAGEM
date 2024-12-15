% Carregar os dados
load('dados.mat');

% Extrair y e u das colunas apropriadas
u = Z(:, 2); % Segunda coluna de Z é o sinal de entrada
y = Z(:, 3); % Terceira coluna de Z é o sinal de saída

% Centralizar os dados
um = u - mean(u); % Centralizar a entrada
ym = y - mean(y); % Centralizar a saída

% Construção da matriz Psim
Psim = [ym(1:end-3), ym(2:end-2), ym(3:end-1), um(2:end-2), um(3:end-1)];

% Cálculo de Ruu e ruy
Ruu = zeros(5, 5);
ruy = zeros(5, 1);

for i = 1:5
    ruy(i) = Psim(:, i)' * ym(4:end);
    for j = 1:5
        Ruu(i, j) = Psim(:, i)' * Psim(:, j);
    end
end

% Estimativa dos parâmetros (tetaChapeu)
tetaChapeu = Ruu \ ruy;

% Estimação de yChapeu (dados de estimação)
yChapeu = y;
yChapeu(4:end) = Psim * tetaChapeu + mean(y);

% Cálculo de resíduos e função de custo
csiMQFC = y - yChapeu;
JMQFC = csiMQFC' * csiMQFC;

% Estimação de yChapeu1 (infinitos passos à frente)
yChapeu1 = ym; % Trabalhar com os dados centralizados
for i = 4:length(yChapeu1)
    yChapeu1(i) = ...
        tetaChapeu(1) * yChapeu1(i-3) + ...
        tetaChapeu(2) * yChapeu1(i-2) + ...
        tetaChapeu(3) * yChapeu1(i-1) + ...
        tetaChapeu(4) * um(i-2) + ...
        tetaChapeu(5) * um(i-1);
end

% Adicionar a média para desfazer a centralização
yChapeu1 = yChapeu1 + mean(y);

% Cálculo de resíduos e função de custo para infinitos passos à frente
csiMQFC1 = y - yChapeu1;
JMQFC1 = csiMQFC1' * csiMQFC1;

% Plot dos resultados
figure;
plot(0:length(y)-1, y, 'k.-', 'LineWidth', 3, 'MarkerSize', 10); hold on; grid on;
plot(3:length(y)-1, yChapeu(4:end), 'r--', 'LineWidth', 3);
plot(0:length(y)-1, yChapeu1, 'g-', 'LineWidth', 3);

% Configuração do gráfico
xlabel('k', 'FontSize', 16);
ylabel('y(k)', 'FontSize', 16);
legend({'Dados reais', 'Um passo à frente (MQ FC)', 'Inf. passos à frente (MQ FC)'}, ...
       'FontSize', 14, 'Location', 'SouthEast');
title(sprintf('MQ com FC: JMQFC = %.2f, JMQFC1 = %.2f', JMQFC, JMQFC1), 'FontSize', 16);
