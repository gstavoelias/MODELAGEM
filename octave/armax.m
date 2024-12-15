% Carregar os dados
clear all;
load('dados.mat');

% Extrair tempo, entrada (u) e saída (y)
t = Z(:, 1); % Coluna 0 é o tempo
u = Z(:, 2); % Coluna 1 é a entrada
y = Z(:, 3); % Coluna 2 é a saída

% Configurações do modelo ARMAX
max_na = 5;
max_nb = 5;
nc = 1; % Ordem da parte de ruído
nk = 1; % Atraso inicial

% Inicializar matrizes de erro e função de custo
eamx = zeros(max_na, max_nb, length(y));
JMQamx = zeros(max_na, max_nb);

% Variáveis para armazenar o melhor modelo
best_na = 0;
best_nb = 0;
best_theta = [];
best_yamx = [];
best_JMQ = Inf;

% Função auxiliar para calcular o modelo ARMAX
function [yamx, theta] = armax_estimation(na, nb, nc, nk, u, y)
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
        % Termos do ruído (opcional, pode ser ajustado)
        for j = 1:nc
            row = [row, 0]; % Placeholder para ruído
        end
        Phi = [Phi; row];
    end

    % Saída correspondente
    y_target = y(max_delay+1:end);

    % Estimativa dos parâmetros usando mínimos quadrados
    theta = pinv(Phi) * y_target;

    % Gerar y estimado
    yamx = zeros(size(y));
    yamx(1:max_delay) = y(1:max_delay);
    for i = max_delay+1:length(y)
        yamx(i) = Phi(i-max_delay, :) * theta;
    end
end

% Loop para varrer valores de na e nb
for na = 1:max_na
    for nb = 1:na % nb limitado por na
        [yamx, theta] = armax_estimation(na, nb, nc, nk, u, y);
        JMQ = sum((y - yamx).^2);
        JMQamx(na, nb) = JMQ;

        % Verificar se este é o melhor modelo
        if JMQ < best_JMQ
            best_na = na;
            best_nb = nb;
            best_theta = theta;
            best_yamx = yamx;
            best_JMQ = JMQ;
        end

        % Plot intermediário
        figure('Position', [100, 100, 1000, 500]);
        plot(t, y, 'k.-', 'LineWidth', 2, 'MarkerSize', 6); hold on;
        plot(t, yamx, 'r--', 'LineWidth', 2);
        grid on;
        xlabel('Tempo', 'FontSize', 14);
        ylabel('y(k)', 'FontSize', 14);
        title(sprintf('Modelo ARMAX (na=%d, nb=%d, JMQ=%.2f)', na, nb, JMQ), 'FontSize', 16);
        legend({'Dados reais', 'ARMAX'}, 'FontSize', 12, 'Location', 'SouthEast');
        hold off;
    end
end

% Plot do melhor modelo
figure('Position', [100, 100, 1000, 500]);
plot(t, y, 'k.-', 'LineWidth', 2, 'MarkerSize', 6); hold on;
plot(t, best_yamx, 'r--', 'LineWidth', 2);
grid on;
xlabel('Tempo', 'FontSize', 14);
ylabel('y(k)', 'FontSize', 14);
title(sprintf('Melhor Modelo ARMAX (na=%d, nb=%d, JMQ=%.2f)', best_na, best_nb, best_JMQ), 'FontSize', 16);
legend({'Dados reais', 'ARMAX'}, 'FontSize', 12, 'Location', 'SouthEast');
hold off;
