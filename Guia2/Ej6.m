% Limpia el espacio de trabajo y cierra todas las figuras abiertas
clear all;
close all;

% Definición de parámetros de la señal
Fs = 1000;                  % Frecuencia de muestreo (Hz)
Ts = 1/Fs;                  % Período de muestreo (s)
t = 0:Ts:1-Ts;              % Vector de tiempo de 1 segundo de duración
N = length(t);              % Número de muestras
k = 400;                    % Parámetro que controla el ancho de la ventana gaussiana
f1 = 50;                    % Frecuencia de la primera señal (Hz)
f2 = 100;                   % Frecuencia de la segunda señal (Hz)

% Creación de la señal compuesta por dos cosenos con diferentes frecuencias
x = zeros(1, N);                % Inicializa la señal x con ceros
half_point = round(N / 2);      % Punto medio para dividir la señal en dos secciones

x(1:half_point) = cos(2*pi*f1*t(1:half_point));    % Primera mitad: coseno de 50 Hz
x(half_point+1:end) = cos(2*pi*f2*t(half_point+1:end)); % Segunda mitad: coseno de 100 Hz

% Inicialización de la matriz para la STFT (Short-Time Fourier Transform)
STFT = zeros(round(N/2), N);       % Matriz para almacenar la STFT

% Vector de frecuencias para la transformación de Fourier
f = 0:1:N-1;

% Cálculo de la STFT utilizando una ventana gaussiana
for n = 1:N
    g = exp(-k*(t - t(n)).^2);               % Ventana gaussiana centrada en t(n)
    x_v = x .* g;                            % Señal modulada por la ventana
    x_v_f = fft(x_v) .* exp(1i*2*pi*t(n).*f); % Transformada de Fourier de la señal modulada
    x_v_f = x_v_f(N/2+1:end);               % Se conserva la mitad positiva del espectro
    STFT(:, n) = x_v_f;                     % Se almacena el resultado en la matriz STFT
end

% Vector de frecuencias para la visualización del espectrograma
k = 0:Fs/N:Fs/2 - Fs/N;

% Visualización de la señal en el dominio del tiempo y del espectrograma
figure;

subplot(2,1,1);                  % Primer gráfico: señal en el dominio del tiempo
plot(t, real(x));
title('Señal en el tiempo');
xlabel('Tiempo (s)');
ylabel('Amplitud');

subplot(2,1,2);                  % Segundo gráfico: espectrograma (STFT)
imagesc(t, k, abs(STFT));         % Visualización de la magnitud de la STFT
set(gca, 'YDir', 'reverse');      % Invierte el eje Y para una mejor interpretación

xlabel('Tiempo (s)');
ylabel('Frecuencia (Hz)');
title('|STFT| ');
colorbar;                         % Barra de color para indicar la magnitud de la STFT
