% Limpia el espacio de trabajo y cierra todas las figuras abiertas
clear all;
close all;

% Definición de parámetros de muestreo
Fs = 1000;               % Frecuencia de muestreo (Hz)
Ts = 1/Fs;               % Período de muestreo (s)
t = 0: Ts :1-Ts;         % Vector de tiempo de 1 segundo de duración
N = length(t);           % Número de muestras
f = -Fs/2: Fs/N: Fs/2 -1; % Vector de frecuencias para la FFT

% Definición de la primera señal gaussiana centrada en t = 0.5
x1 = exp(-2000*(t-0.5).^2);

% Representación en el dominio del tiempo de x1
figure;
subplot(321);
plot(t, x1);
xlabel('Tiempo(s)');
ylabel('Amplitud');
title('x1(t)');

% Cálculo y visualización de la FFT de x1
fft_x1 = fft(x1);           % Transformada de Fourier
fft_x1 = fftshift(fft_x1);  % Centra el cero en la frecuencia
fft_x1 = abs(fft_x1);       % Módulo de la FFT

subplot(322);
plot(f, fft_x1);
title('|x1(f)|');
xlabel('Frecuencia (Hz)');
ylabel('Amplitud');

% Definición de la segunda señal gaussiana centrada en t = 0.25
x2 = exp(-2000*(t-0.25).^2);

% Representación en el dominio del tiempo de x2
subplot(323);
plot(t, x2);
title('x2(t)');
xlabel('Tiempo(s)');
ylabel('Amplitud');

% Cálculo y visualización de la FFT de x2
fft_x2 = fft(x2);
fft_x2 = fftshift(fft_x2);
fft_x2 = abs(fft_x2);

subplot(324);
plot(f, fft_x2);
title('|x2(f)|');
xlabel('Frecuencia (Hz)');
ylabel('Amplitud');

% Definición de la tercera señal gaussiana centrada en t = 0.75
x3 = exp(-2000*(t-0.75).^2);

% Representación en el dominio del tiempo de x3
subplot(325);
plot(t, x3);
title('x3(t)');
xlabel('Tiempo(s)');
ylabel('Amplitud');

% Cálculo y visualización de la FFT de x3
fft_x3 = fft(x3);
fft_x3 = fftshift(fft_x3);
fft_x3 = abs(fft_x3);

subplot(326);
plot(f, fft_x3);
title('|x3(f)|');
xlabel('Frecuencia (Hz)');
ylabel('Amplitud');


% Comentario:
% Al desplazar las señales temporalmente, se observa que el módulo de su
% transformada de Fourier se mantiene invariante, lo que refleja la propiedad
% de la invariancia del módulo de la FFT frente a desplazamientos temporales.


