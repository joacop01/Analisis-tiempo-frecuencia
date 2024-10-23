Fs = 1000;                  % Frecuencia de muestreo
Ts = 1/Fs;                  % Intervalo de muestreo
t = 0:Ts:1-Ts;              % Vector de tiempo
N = length(t);              % Longitud de la señal
k = 200;                    % Parámetro de la ventana gaussiana
len = 0.1:0.01:1;           % Diferentes longitudes de la ventana
sigma = 18 ./ len.^2;       % Parámetro sigma para cada longitud de ventana
f1 = 100;                   % Frecuencia de la primera señal
f2 = 150;                   % Frecuencia de la segunda señal
alpha1 = 2;                 % Parámetro alpha para la entropía de Rényi (caso 1)
alpha2 = 3;                 % Parámetro alpha para la entropía de Rényi (caso 2)

% x1 = cos(2*pi*f1*t + 2*pi*100*t.^2);  % Primera señal
% x2 = cos(2*pi*f2*t + 2*pi*100*t.^2);  % Segunda señal
x1 = cos(2*pi*f1*t);                  % Primera señal
x2 = cos(2*pi*f2*t);                  % Segunda señal
x = x1 + x2;                          % Señal combinada

P_signal = mean(x.^2);                % Potencia de la señal original

% Generar ruido gaussiano con la misma potencia que la señal
ruido_gaussiano = sqrt(P_signal) * randn(size(x));  % Escalar el ruido para que tenga la misma potencia
x_ruido = x + ruido_gaussiano;               % Señal con ruido

% Parámetros para la STFT
STFT = zeros(N, N);                    % Matriz para almacenar el resultado de la STFT

% Inicialización de los vectores para almacenar la entropía de Rényi
renyi_alpha1 = zeros(1, length(len)); 
renyi_alpha2 = zeros(1, length(len)); 

for l = 1:length(len)
    for n = 1:N
        a = n / Fs;                           % Desplazamiento de la ventana
        g = exp(-sigma(l) * (t - a).^2);      % Ventana gaussiana
        x_v = x_ruido .* g;                         % Multiplicación de la señal por la ventana
        x_v_f = fftshift(fft(x_v));           % FFT de la señal segmentada
        STFT(:, n) = x_v_f;                   % Almacenar en la matriz de STFT
    end

    % Cálculo del espectrograma
    spec = abs(STFT).^2;                      % Potencia espectral (magnitud al cuadrado)

    % Normalización del espectrograma
    spec_norm = spec ./ sum(sum(spec));         % Normalización para que el espectrograma sea una distribución de probabilidad

    % Cálculo de la entropía de Rényi para alpha = 2
    renyi_alpha1(l) = (1 / (1 - alpha1)) * log2(sum(sum(spec_norm.^alpha1)));

    % Cálculo de la entropía de Rényi para alpha = 3
    renyi_alpha2(l) = (1 / (1 - alpha2)) * log2(sum(sum(spec_norm.^alpha2)));

end

subplot(2, 1, 1);
plot(t, x, 'LineWidth', 1.5);
title('Señal Original');
xlabel('Tiempo (s)');
ylabel('Amplitud');
grid on;

% Subplot para la señal ruidosa con SNR = 0 dB
subplot(2, 1, 2);
plot(t, x_ruido, 'LineWidth', 1.5);
title('Señal Ruidosa (SNR = 0 dB)');
xlabel('Tiempo (s)');
ylabel('Amplitud');
grid on;


% Graficar la entropía de Rényi en función de la longitud de la ventana
figure;
plot(len, renyi_alpha1, '-o', 'DisplayName', ['α = ', num2str(alpha1)]);
hold on;
plot(len, renyi_alpha2, '-s', 'DisplayName', ['α = ', num2str(alpha2)]);
xlabel('Longitud de la ventana');
ylabel('Entropía de Rényi');
title('Entropía de Rényi vs Longitud de la ventana');
legend show;

% Encontrar el mínimo global para alpha = 2 y alpha = 3
[min_renyi1, index_min_renyi1] = min(renyi_alpha1);
[min_renyi2, index_min_renyi2] = min(renyi_alpha2);

% Mostrar resultados del mínimo
disp(['Mínimo global de la entropía de Rényi (α = ', num2str(alpha1), '): ', num2str(min_renyi1), ' en longitud de ventana: ', num2str(len(index_min_renyi1)), 's']);
disp(['Mínimo global de la entropía de Rényi (α = ', num2str(alpha2), '): ', num2str(min_renyi2), ' en longitud de ventana: ', num2str(len(index_min_renyi2)),'s']);