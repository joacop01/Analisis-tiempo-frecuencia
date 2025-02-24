% Limpiar el espacio de trabajo y cerrar figuras abiertas
clear all;
close all;

% Definición de parámetros de muestreo
Fs = 1000;                 % Frecuencia de muestreo (Hz)
Ts = 1/Fs;                 % Periodo de muestreo (s)
t = 0:Ts:1-Ts;             % Vector de tiempo de 1 segundo de duración
N = length(t);             % Número de muestras
k = (-Fs/2:Fs/N:Fs/2 - Fs/N);  % Vector de frecuencias centrado en cero

% Frecuencias de las señales
f1 = 200;                  % Frecuencia de la primera señal (Hz)
f2 = 100;                  % Frecuencia de la segunda señal (Hz)

% Inicialización de la señal
x = zeros(1, N);           % Vector de ceros para la señal
half_point = round(N / 2); % Punto medio del vector de tiempo

% Construcción de la señal con dos frecuencias diferentes en la primera y segunda mitad
x(1:half_point) = cos(2*pi*f1*t(1:half_point));        % Primera mitad con frecuencia f1
x(half_point+1:end) = cos(2*pi*f2*t(half_point+1:end)); % Segunda mitad con frecuencia f2

% Transformada de Fourier de la señal y centrado del espectro
x_f = fftshift(fft(x));    % Aplicación de la FFT y desplazamiento del espectro

% Inicialización de la matriz para la STFT (Short-Time Fourier Transform)
STFT = zeros(N, N);        % Matriz para almacenar los resultados de la STFT
l = length(t) / 2;         % Mitad de la longitud de la señal
T = 0.2;                   % Ancho de la ventana temporal

% Cálculo de la STFT utilizando una ventana basada en funciones sinc
for f = 1:N
    % Definición de la ventana espectral usando combinaciones de funciones sinc
    g = 1/2 * T * sinc((k-k(f))*T) + 1/4*T * sinc(((k-k(f))-1/T) *T) + 1/4*T * sinc((k-k(f)+1/T) *T);
    
    % Filtrado de la señal en el dominio de la frecuencia
    x_v = x_f .* g;                           
    
    % Transformada inversa de Fourier para obtener la señal filtrada en el dominio temporal
    STFT(f, :) = ifft(x_v);        
end

% Definición del vector de frecuencias positivas
k_2 = 0: Fs/N: Fs/2 -1;

% Visualización de los resultados
figure;

% Representación de la señal en el dominio de la frecuencia
subplot(2,1,1);
plot(k_2, abs(x_f(N/2+1:end)));
title('Señal en la frecuencia');
xlabel('Frecuencia (Hz)');
ylabel('Amplitud');

% Representación del espectrograma (STFT)
subplot(2,1,2);
imagesc(t, k_2, abs(STFT(N/2+1:end,:)));  % Magnitud de la STFT

% Invertir el eje de frecuencias para mostrar frecuencias más altas abajo
set(gca, 'YDir', 'reverse'); 

xlabel('Tiempo (s)');
ylabel('Frecuencia (Hz)');
title('|STFT|');
colorbar;
