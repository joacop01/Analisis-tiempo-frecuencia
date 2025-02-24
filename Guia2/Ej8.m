clear all;
close all;

% Parámetros de la señal
Fs = 1000;                  % Frecuencia de muestreo (Hz)
Ts = 1/Fs;                  % Período de muestreo (s)
t = 0:Ts:1-Ts;              % Vector de tiempo
N = length(t);              % Número de muestras
T = 0.1;                    % Ancho de la ventana (s)

f1 = 50;                    % Frecuencia de la primera señal (Hz)
f2 = 100;                   % Frecuencia de la segunda señal (Hz)

% Generación de la señal
x = zeros(1, N);            % Inicialización de la señal
half_point = round(N / 2);  % Punto medio de la señal

% Combinación de dos cosenos en diferentes mitades del tiempo
x(1:half_point) = cos(2*pi*f1*t(1:half_point));
x(half_point+1:end) = cos(2*pi*f2*t(half_point+1:end));

% Inicialización de la STFT
STFT = zeros(N/2, N);

f = (-Fs/2:Fs/N:Fs/2 - Fs/N);  % Vector de frecuencias centradas

% Cálculo de la STFT con una ventana de coseno
for n = 1:N
    chi = (t - t(n) >= -T/2) & (t - t(n) <= T/2);  % Ventana rectangular
    g = (1 + cos(2*pi*(t-t(n))/T)) .* chi;         

    x_v = x .* g;                                 % Señal ventaneada
    x_v_f = fft(x_v) .* exp(1i*2*pi*t(n).*f);     % Transformada de Fourier
    x_v_f = x_v_f(N/2+1:end);                     % Solo la mitad positiva
    STFT(:, n) = x_v_f;                           % Guardar en la matriz STFT
end

k = 0:Fs/N:Fs/2 - Fs/N;  % Vector de frecuencias positivas

% Gráficos
figure;

subplot(2,1,1);
plot(t, x);
title('Señal en el tiempo');
xlabel('Tiempo (s)');
ylabel('Amplitud');

subplot(2,1,2);
imagesc(t, k, abs(STFT));   % Espectrograma
set(gca, 'YDir', 'normal'); % Mantener el eje de frecuencia ascendente

xlabel('Tiempo (s)');
ylabel('Frecuencia (Hz)');
title('|STFT|');
colorbar;
