clear all; % Borra todas las variables del workspace
close all; % Cierra todas las figuras abiertas

Fs = 1000; % Frecuencia de muestreo en Hz
Ts = 1/Fs; % Período de muestreo en segundos
t = 0:Ts:1-Ts; % Vector de tiempo de 1 segundo de duración con paso Ts
N = length(t); % Número total de muestras
f = -Fs/2: Fs/N: Fs/2 -1; % Vector de frecuencias centrado en 0 Hz

% Generación de una señal gaussiana en el dominio del tiempo
x1 = exp(-2000*(t-0.5).^2); % Pulso gaussiano centrado en t = 0.5s

% Cálculo de la Transformada de Fourier de x1
fft_x1 = fft(x1); % Cálculo de la FFT
fft_x1 = fftshift(fft_x1); % Desplazamiento del espectro para centrarlo en 0 Hz
fft_x1 = abs(fft_x1); % Se toma el valor absoluto para obtener la magnitud

% Modulación de x1 con una onda coseno de 250 Hz
x4 = x1 .* cos(2*pi*250*t); % Multiplicación en el dominio del tiempo

figure; % Nueva figura para graficar

% Gráfico de la señal original x1 en el dominio del tiempo
subplot(321);
plot(t, x1);
title('x1(t)');
xlabel('Tiempo (s)');
ylabel('Amplitud');

% Gráfico del espectro de x1 en el dominio de la frecuencia
subplot(322);
plot(f, fft_x1);
title('|x1(f)|');
xlabel('Frecuencia (Hz)');
ylabel('Amplitud');

% Gráfico de la señal modulada x4 en el dominio del tiempo
subplot(323);
plot(t, x4);
title('x4(t)');
xlabel('Tiempo (s)');
ylabel('Amplitud');

% Cálculo de la FFT de x4
fft_x4 = fft(x4);
fft_x4 = fftshift(fft_x4);
fft_x4 = abs(fft_x4);

% Gráfico del espectro de x4 en el dominio de la frecuencia
subplot(324);
plot(f, fft_x4);
title('|x4(f)|');
xlabel('Frecuencia (Hz)');
ylabel('Amplitud');

% Modulación de x1 con una exponencial compleja de 350 Hz
x5 = x1 .* exp(1i*2*pi*350*t); % Modulación en frecuencia

% Gráfico de la señal modulada x5 en el dominio del tiempo
subplot(325);
hold on;
plot(t, real(x5), 'r'); % Parte real en rojo
title('x5(t)');
plot(t, imag(x5), 'b'); % Parte imaginaria en azul
legend('Re(x5(t))', 'Im(x5(t))');
xlabel('Tiempo (s)');
ylabel('Amplitud');
hold off;

% Cálculo de la FFT de x5
fft_x5 = fft(x5);
fft_x5 = fftshift(fft_x5);
fft_x5 = abs(fft_x5);

% Gráfico del espectro de x5 en el dominio de la frecuencia
subplot(326);
plot(f, fft_x5);
title('|x5(f)|');
xlabel('Frecuencia (Hz)');
ylabel('Amplitud');

% Comentarios finales:
% - La multiplicación por un coseno desplaza el espectro de la señal en ± la frecuencia del coseno,
%   lo que se observa en el gráfico de |x4(f)|.
% - La multiplicación por una exponencial compleja genera un desplazamiento en una sola dirección,
%   eliminando las frecuencias negativas, lo que se observa en |x5(f)|.