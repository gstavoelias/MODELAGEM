% Carregar o arquivo .mat
load('dados.mat');

% Definir as variáveis u e y a partir das colunas de Z
u = Z(:, 2); % Segunda coluna de Z é o sinal de entrada
y = Z(:, 3); % Terceira coluna de Z é o sinal de saída

% Verificar tamanho mínimo de y e u
if length(y) < 4 || length(u) < 3
    error('Os vetores y e u precisam ter pelo menos 4 e 3 elementos, respectivamente.');
end

% Construir a matriz Psi e o vetor y1
Psi = [y(1:end-3), y(2:end-2), y(3:end-1), u(2:end-2), u(3:end-1)];
y1 = y(4:end);

% Estimar os parâmetros (tetaChapeu)
tetaChapeu = (Psi' * Psi) \ (Psi' * y1);

% Estimação de yChapeu (um passo à frente)
yChapeu = y;
yChapeu(4:end) = Psi * tetaChapeu;
csiMQ = y - yChapeu;
JMQ = csiMQ' * csiMQ;

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
csiMQ1 = y - yChapeu1;
JMQ1 = csiMQ1' * csiMQ1;

% Plot dos resultados
figure;
plot(0:length(y)-1, y, 'k.-', 'LineWidth', 3, 'MarkerSize', 10); hold on; grid on;
plot(0:length(yChapeu)-1, yChapeu, 'r.-', 'LineWidth', 3, 'MarkerSize', 10);
plot(0:length(yChapeu1)-1, yChapeu1, 'g-', 'LineWidth', 3);

% Configuração do gráfico
xlabel('k', 'FontSize', 16);
ylabel('y(k)', 'FontSize', 16);
legend({'Dados de estimação', 'Um passo à frente', 'Inf. passos à frente'}, ...
       'FontSize', 14, 'Location', 'SouthEast');
title(sprintf('MQP: JMQ = %.2f, JMQ1 = %.2f', JMQ, JMQ1), 'FontSize', 16);
