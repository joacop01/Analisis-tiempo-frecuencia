Fs = 1000;                  % Frecuencia de muestreo
Ts = 1/Fs;                  % Intervalo de muestreo
t = 0:Ts:1-Ts;              % Vector de tiempo
N = length(t);              % Longitud de la señal
% Vector de frecuencias
k = (-Fs/2:Fs/N:Fs/2 - Fs/N);               % Vector de frecuencias centrado

sigma = 300;                   % Parámetro de la ventana gaussiana
f1 = 50;                    % Frecuencia de la señal
f2 = 100;                   % Frecuencia de la señal
% Crear la señal

x = zeros(1, N);                % Inicializar la señal con ceros
half_point = round(N / 2);       % Punto intermedio del tiempo

% Primera mitad de la señal (10 Hz)
x(1:half_point) = cos(2*pi*f1*t(1:half_point));

% Segunda mitad de la señal (30 Hz)
x(half_point+1:end) = cos(2*pi*f2*t(half_point+1:end));

x_f = fftshift(fft(x));
% Parámetros para la STFT
STFT = zeros(N, N);         % Matriz para almacenar el resultado de la STFT
l = length(t) / 2;          % Ventana centrada


% Calcular la STFT
for f = 1:N
    % Desplazamiento de la ventana
    g = sqrt(pi/sigma)*exp((-pi*(k-k(f)).^2)/sigma);                 % Ventana gaussiana
    x_v = x_f .* g;                           % Multiplicación de la señal por la ventana
    STFT(f, :) = ifft(x_v);        % FFT y fftshift para centrar la frecuencia
end

k_2 = 0: Fs/N: Fs/2 -1;
%Graficar la señal original y su STFT
figure;
subplot(2,1,1);
plot(k, abs(x_v));
hold on;
plot(k, g);
title('Señal en la frecuencia');
xlabel('Frecuencia (Hz)');
ylabel('Amplitud');

% Mostrar el espectrograma (STFT en función del tiempo y frecuencia)
subplot(2,1,2);
imagesc(t, k_2, abs(STFT(N/2+1:end, :)));   % Magnitud de la STFT

% Invertir el eje de frecuencias (frecuencias más altas abajo)
set(gca, 'YDir', 'reverse');  % Ahora las frecuencias altas estarán abajo

xlabel('Tiempo (s)');
ylabel('Frecuencia (Hz)');
title('|STFT| ');
colorbar;