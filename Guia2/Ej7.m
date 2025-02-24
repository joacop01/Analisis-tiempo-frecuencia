clear all;          % Limpia todas las variables del espacio de trabajo
close all;          % Cierra todas las figuras abiertas

Fs = 1000;          % Frecuencia de muestreo en Hz
Ts = 1/Fs;          % Período de muestreo en segundos
t = 0:Ts:1-Ts;      % Vector de tiempo de 1 segundo de duración
N = length(t);      % Número total de muestras

k = (-Fs/2:Fs/N:Fs/2 - Fs/N);  % Vector de frecuencias centrado en cero

sigma = 300;        % Parámetro de dispersión para la ventana gaussiana
f1 = 50;            % Frecuencia de la primera señal (en Hz)
f2 = 100;           % Frecuencia de la segunda señal (en Hz)

x = zeros(1, N);    % Inicialización de la señal
half_point = round(N / 2);  % Punto medio del vector de tiempo

% Construcción de la señal: combinación de dos cosenos de diferente frecuencia
x(1:half_point) = cos(2*pi*f1*t(1:half_point));       % Primera mitad con f1
x(half_point+1:end) = cos(2*pi*f2*t(half_point+1:end)); % Segunda mitad con f2

x_f = (fft(x));     % Transformada de Fourier de la señal completa

STFT = zeros(N, N); % Inicialización de la matriz para la STFT

% Cálculo de la STFT mediante una ventana gaussiana
for f = 1:N
    g = sqrt(pi/sigma) * exp((-pi * (k - k(f)).^2) / sigma); % Ventana gaussiana centrada en k(f)
    x_v = x_f .* g;         % Aplicación de la ventana en el dominio de la frecuencia
    STFT(f, :) = ifft(x_v); % Transformada inversa para obtener la STFT en el tiempo
end

k_2 = 0: Fs/N: Fs/2 - 1;    % Vector de frecuencias positivas para la visualización

% Gráficos de resultados
figure;
subplot(2,1,1);             % Primer subplot: dominio de la frecuencia
plot(k(N/2:end), x_f(N/2:end)); % Espectro de la señal (parte positiva)
title('Señal en la frecuencia');
xlabel('Frecuencia (Hz)');
ylabel('Amplitud');

subplot(2,1,2);             % Segundo subplot: espectrograma (STFT)
imagesc(t, k_2, abs(STFT(N/2+1:end, :))); % Representación de la magnitud de la STFT
set(gca, 'YDir', 'reverse'); % Invertir el eje Y para mostrar frecuencias más bajas arriba

xlabel('Tiempo (s)');       % Etiqueta del eje X
ylabel('Frecuencia (Hz)');  % Etiqueta del eje Y
title('|STFT|');            % Título del gráfico
colorbar;                   % Barra de colores para indicar la magnitud
