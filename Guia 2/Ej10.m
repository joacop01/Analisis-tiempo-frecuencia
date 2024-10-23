addpath('C:\Repositorios\Curso Análisis Tiempo-Frecuencia\Funciones');
Fs = 1000;
Ts = 1/Fs;
t = 0: Ts: 1-Ts;

x1 = cos(2*pi*100*t + 2*pi*100*t.^2);
x2 = cos(2*pi*150*t + 2*pi*100*t.^2);
x = x1 + x2;

% Calcular la potencia de la señal
P_signal = mean(x.^2);  % Potencia de la señal original
% Generar ruido gaussiano con la misma potencia que la señal
ruido_gaussiano = sqrt(P_signal) * randn(size(x));  % Escalar el ruido para que tenga la misma potencia
P_noise = mean(ruido_gaussiano.^2);  % Potencia del ruido
SNR = 10 * log10(P_signal / P_noise);
% Mostrar el valor de la SNR (debería ser cercano a 0 dB)
disp(['Relación señal-ruido (SNR): ', num2str(SNR), ' dB']);

figure;
subplot(2,1,1);
plot(t, x);
title('Señal Original');
xlabel('Tiempo (s)');
ylabel('Amplitud');

x = x +ruido_gaussiano;

subplot(2,1,2);
plot(t, x);
title('Señal ruidosa');
xlabel('Tiempo (s)');
ylabel('Amplitud');

STFT_gauss = STFT_Gauss(x, t, 200);
STFT_hann = STFT_Hann(x, t, 0.2);

k = 0:Fs/N: Fs/2 - Fs/N;                    % Frecuencias positivas

% Mostrar el espectrograma (STFT en función del tiempo y frecuencia)
figure;
subplot(211);
imagesc(t, k, abs(STFT_gauss));   % Magnitud de la STFT
% Invertir el eje de frecuencias (frecuencias más altas abajo)
set(gca, 'YDir', 'reverse');  % Ahora las frecuencias altas estarán abajo
xlabel('Tiempo (s)');
ylabel('Frecuencia (Hz)');
title('Espectrograma STFT con ventana Gaussiana');
colorbar;

% Mostrar el espectrograma (STFT en función del tiempo y frecuencia)
subplot(212);
imagesc(t, k, abs(STFT_hann));   % Magnitud de la STFT
% Invertir el eje de frecuencias (frecuencias más altas abajo)
set(gca, 'YDir', 'reverse');  % Ahora las frecuencias altas estarán abajo
xlabel('Tiempo (s)');
ylabel('Frecuencia (Hz)');
title('Espectrograma STFT con ventana de Hann');
colorbar;
