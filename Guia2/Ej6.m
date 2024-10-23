Fs = 1000;                  % Frecuencia de muestreo
Ts = 1/Fs;                  % Intervalo de muestreo
t = 0:Ts:1-Ts;              % Vector de tiempo
N = length(t);              % Longitud de la señal
k = 18;                   % Parámetro de la ventana gaussiana
f1 = 50;                    % Frecuencia de la señal
f2 = 100;                   % Frecuencia de la señal
f = 100;
% x = exp(1i*2*pi*f*t);        % Señal compleja
x = zeros(1, N);                % Inicializar la señal con ceros
half_point = round(N / 2);       % Punto intermedio del tiempo

% Primera mitad de la señal (10 Hz)
x(1:half_point) = cos(2*pi*f1*t(1:half_point));

% Segunda mitad de la señal (30 Hz)
x(half_point+1:end) = cos(2*pi*f2*t(half_point+1:end));

% Parámetros para la STFT
STFT = zeros(N/2, N);        % Matriz para almacenar el resultado de la STFT

% Vector de frecuencias
f = (-Fs/2:Fs/N:Fs/2 - Fs/N); % Vector de frecuencias centrado

% Calcular la STFT
for n = 1:N
    a = n/ Fs;                           % Desplazamiento de la ventana
    g = exp(-k*(t - a).^2);                 % Ventana gaussiana
    x_v = x .* g;                           % Multiplicación de la señal por la ventana
    x_v_f = fftshift(fft(x_v)).*exp(1i*2*pi*t/N);
    x_v_f = x_v_f(N/2+1:end);
    STFT(:, n) = x_v_f;                     % FFT y fftshift para centrar la frecuencia
end

k = 0:Fs/N: Fs/2 - Fs/N;                    % Frecuencias positivas

% Graficar la señal original y su STFT
figure;
subplot(2,1,1);
plot(t, real(x));
title('Señal en el tiempo');
xlabel('Tiempo (s)');
ylabel('Amplitud');
%hold on;
%plot(t, g);
%title('Ventana gaussiana');
%hold off;

% Mostrar el espectrograma (STFT en función del tiempo y frecuencia)
subplot(2,1,2);
imagesc(t, k, abs(STFT));   % Magnitud de la STFT

% Invertir el eje de frecuencias (frecuencias más altas abajo)
set(gca, 'YDir', 'reverse');  % Ahora las frecuencias altas estarán abajo

xlabel('Tiempo (s)');
ylabel('Frecuencia (Hz)');
title('|STFT| ');
colorbar;