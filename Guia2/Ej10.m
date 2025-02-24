% Agregar la carpeta de funciones al path de MATLAB
addpath('C:\Repositorios\Analisis Tiempo Frecuencia y Descomposicion de Señales\Funciones');

% Definición de parámetros de muestreo
Fs = 1000;         % Frecuencia de muestreo (Hz)
Ts = 1/Fs;         % Período de muestreo (s)
t = 0:Ts:1-Ts;     % Vector de tiempo de 1 segundo

% Generación de dos señales cosenoidales con modulación en frecuencia
x1 = cos(2*pi*100*t + 2*pi*100*t.^2);  % Señal 1: frecuencia base 100 Hz
x2 = cos(2*pi*150*t + 2*pi*100*t.^2);  % Señal 2: frecuencia base 150 Hz
x = x1 + x2;                           % Señal compuesta

% Calcular la potencia de la señal original
P_signal = mean(x.^2);  % Potencia de la señal (media del cuadrado de la amplitud)

% Generar ruido gaussiano con la misma potencia que la señal
ruido_gaussiano = sqrt(P_signal) * randn(size(x));  % Ruido blanco gaussiano
P_noise = mean(ruido_gaussiano.^2);                % Potencia del ruido

% Calcular la relación señal-ruido (SNR)
SNR = 10 * log10(P_signal / P_noise);  % SNR en decibelios (dB)
disp(['Relación señal-ruido (SNR): ', num2str(SNR), ' dB']);

% Visualización de la señal original y la señal ruidosa
figure;
subplot(2,1,1);
plot(t, x);
title('Señal Original');
xlabel('Tiempo (s)');
ylabel('Amplitud');

% Añadir ruido gaussiano a la señal original
x = x + ruido_gaussiano;

subplot(2,1,2);
plot(t, x);
title('Señal ruidosa');
xlabel('Tiempo (s)');
ylabel('Amplitud');

% Cálculo de la Transformada de Fourier de Tiempo Corto (STFT)
STFT_gauss = STFT_Gauss(x, t, 200);  % STFT usando ventana Gaussiana 
STFT_hann = STFT_Hann(x, t, 0.2);    % STFT usando ventana de Hann

% Definición del vector de frecuencias positivas
N = length(t);                       % Número de muestras
k = 0:Fs/N:Fs/2 - Fs/N;              % Vector de frecuencias (Hz)

% Visualización de los espectrogramas (magnitud de la STFT)
figure;

% Espectrograma con ventana Gaussiana
subplot(2,1,1);
imagesc(t, k, abs(STFT_gauss));      % Magnitud de la STFT
set(gca, 'YDir', 'reverse');         % Invertir eje Y para mostrar frecuencias altas abajo
xlabel('Tiempo (s)');
ylabel('Frecuencia (Hz)');
title('Espectrograma STFT con ventana Gaussiana');
colorbar;                            % Barra de color para la magnitud

% Espectrograma con ventana de Hann
subplot(2,1,2);
imagesc(t, k, abs(STFT_hann));       % Magnitud de la STFT
set(gca, 'YDir', 'reverse');         % Invertir eje Y para mostrar frecuencias altas abajo
xlabel('Tiempo (s)');
ylabel('Frecuencia (Hz)');
title('Espectrograma STFT con ventana de Hann');
colorbar;                            % Barra de color para la magnitud


%Comentario:
%Se observa que con ambas STFT (ventana Gaussiana y ventana de Hann) se
%obtiene un resultado muy similar.