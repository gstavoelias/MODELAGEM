% Carregar os dados
clear all;
load('dados.mat');

% Extrair tempo, entrada (u) e saída (y)
t = Z(:, 1); % Coluna 0 é o tempo
u = Z(:, 2); % Coluna 1 é a entrada
y = Z(:, 3); % Coluna 2 é a saída

% Configurações do modelo ARX
max_na = 5;
max_nb = 5;
nk = 1; % Atraso inicial

% Inicializar matrizes de erro e função de custo
earx = zeros(max_na, max_nb, length(y));
JMQarx = zeros(max_na, max_nb);

% Função auxiliar para calcular o modelo ARX
function [yarx, theta] = arx_estimation(na, nb, nk, u, y)
    max_delay = max(na, nb + nk - 1);
    Phi = [];
    % Construção da matriz de regressores
    for i = max_delay+1:length(y)
        row = [];
        % Termos autoregressivos de y
        for j = 1:na
            row = [row, -y(i-j)];
        end
        % Termos de regressão da entrada u
        for j = 1:nb
            row = [row, u(i-j-nk+1)];
        end
        Phi = [Phi; row];
    end

    % Saída alvo
    y_target = y(max_delay+1:end);

    % Estimativa dos parâmetros usando mínimos quadrados
    theta = pinv(Phi) * y_target;

    % Gerar y estimado
    yarx = zeros(size(y));
    yarx(1:max_delay) = y(1:max_delay);
    for i = max_delay+1:length(y)
        yarx(i) = Phi(i-max_delay, :) * theta;
    end
end

% Loop para varrer valores de na e nb
for na = 1:max_na
    for nb = 1:na % nb limitado por na
        [yarx, ~] = arx_estimation(na, nb, nk, u, y);
        earx(na, nb, :) = y - yarx;
        JMQarx(na, nb) = sum((y - yarx).^2);

        % Plot intermediário
        figure('Position', [100, 100, 1000, 500]);
        plot(t, y, 'k.-', 'LineWidth', 2, 'MarkerSize', 6); hold on;
        plot(t, yarx, 'r--', 'LineWidth', 2);
        grid on;
        xlabel('Tempo', 'FontSize', 14);
        ylabel('y(k)', 'FontSize', 14);
        title(sprintf('Modelo ARX (na=%d, nb=%d, JMQ=%.2f)', na, nb, JMQarx(na, nb)), 'FontSize', 16);
        legend({'Dados reais', 'ARX'}, 'FontSize', 12, 'Location', 'SouthEast');
        hold off;
    end
end

% Simulação com valores fixos de na e nb
na = 5;
nb = 1;
[yarx, theta] = arx_estimation(na, nb, nk, u, y);

% Plot final
figure('Position', [100, 100, 1000, 500]);
plot(t, y, 'k.-', 'LineWidth', 2, 'MarkerSize', 6); hold on;
plot(t, yarx, 'r--', 'LineWidth', 2);
grid on;
xlabel('Tempo', 'FontSize', 14);
ylabel('y(k)', 'FontSize', 14);
title(sprintf('Modelo ARX Final (na=%d, nb=%d, JMQ=%.2f)', na, nb, sum((y - yarx).^2)), 'FontSize', 16);
legend({'Dados reais', 'ARX'}, 'FontSize', 12, 'Location', 'SouthEast');
hold off;
