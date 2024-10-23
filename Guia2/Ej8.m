
% 
% % Condición para que t esté en [-T/2, T/2]

% 
% % La función x(t)
% 
% 
% plot(t, x);

Fs = 1000;                  % Frecuencia de muestreo
Ts = 1/Fs;                  % Intervalo de muestreo
t = 0:Ts:1-Ts;              % Vector de tiempo
N = length(t);              % Longitud de la señal
T = 0.01;
f1 = 50;                    % Frecuencia de la señal
f2 = 100;                   % Frecuencia de la señal
chi = (t >= -T/2) & (t <= T/2);
% Crear la señal
x = zeros(1, N);                % Inicializar la señal con ceros
half_point = round(N / 2);       % Punto intermedio del tiempo

% Primera mitad de la señal (10 Hz)
% x(1:half_point) = cos(2*pi*f1*t(1:half_point));
% 
% % Segunda mitad de la señal (30 Hz)
% x(half_point+1:end) = cos(2*pi*f2*t(half_point+1:end));
x = dirac(t-(N/2)*Ts);
% Parámetros para la STFT
STFT = zeros(N/2, N);         % Matriz para almacenar el resultado de la STFT
l = length(t) / 2;          % Ventana centrada

% Vector de frecuencias
f = (-Fs/2:Fs/N:Fs/2 - Fs/N);               % Vector de frecuencias centrado

% Calcular la STFT
for n = 1:N
% Desplazamiento de la ventana
    chi = (t - n/Fs >= -T/2) & (t - n/Fs <= T/2);
    g = (1 + cos(2*pi*(t-n/Fs)/T)) .* chi;
    x_v = x .* g;                           % Multiplicación de la señal por la ventana
    x_v_f = fftshift(fft(x_v)).*exp(1i*2*pi*t/N);
    x_v_f = x_v_f(N/2+1:end);
    STFT(:, n) = x_v_f;                     % FFT y fftshift para centrar la frecuencia
end

figure;
plot(t, g);
hold on;
plot(t, x_v);
hold off;

k = 0:Fs/N: Fs/2 - Fs/N;                    % Frecuencias positivas
% Graficar la señal original y su STFT
figure;
subplot(2,1,1);
plot(t, x);
title('Señal en el tiempo');
xlabel('Tiempo (s)');
ylabel('Amplitud');

% Mostrar el espectrograma (STFT en función del tiempo y frecuencia)
subplot(2,1,2);
imagesc(t, k, abs(STFT));   % Magnitud de la STFT

% Invertir el eje de frecuencias (frecuencias más altas abajo)
set(gca, 'YDir', 'reverse');  % Ahora las frecuencias altas estarán abajo

xlabel('Tiempo (s)');
ylabel('Frecuencia (Hz)');
title('|STFT| ');
colorbar;