clear all;
close all;
addpath('../Funciones');
Fs = 1000;             % Frecuencia de muestreo
Ts = 1/Fs;             % Periodo de muestreo
t = 0:Ts:1-Ts;         % Vector de tiempo
N = length(t);         % Número de puntos en el tiempo
f = -Fs/2 : Fs/N : Fs/2-Fs/N;  % Vector de frecuencias


x1 = cos(2*pi*100*t + 2*pi*100*t.^2);
x2 = cos(2*pi*150*t + 2*pi*100*t.^2);

x = x1 ;%+ x2;

% Calcular la potencia de la señal
P_signal = mean(x.^2);  % Potencia de la señal original

% Generar ruido gaussiano con la misma potencia que la señal
ruido_gaussiano = sqrt(P_signal) * randn(size(x));  % Escalar el ruido para que tenga la misma potencia

% Graficar la señal original y la señal con ruido

% Calcular la SNR resultante
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

% x = x + ruido_gaussiano;

subplot(2,1,2);
plot(t, x);
title('Señal ruidosa');
xlabel('Tiempo (s)');
ylabel('Amplitud');

scales = 0.0001:0.01:0.5;

% Inicializar la matriz de CWT
CWT_mex = CWT_mexhat(x, f, 1000, scales);
CWT_mor = CWT_morlet(x, f, 800, scales, 10);
CWT_bum = CWT_bump(x, f, 8, scales, 10);

%Graficar la CWT
figure;
subplot(311);
imagesc(t, scales, abs(CWT_mex));  % Graficar la magnitud de la CWT
xlabel('Tiempo');
ylabel('Escalas');
title('CWT con ondita mexhat (analitica)');
colorbar;

%Graficar la CWT
subplot(312);
imagesc(t, scales, abs(CWT_bum));  % Graficar la magnitud de la CWT
xlabel('Tiempo');
ylabel('Escalas');
title('CWT con ondita bump (analitica)');
colorbar;

%Graficar la CWT
subplot(313);
imagesc(t, scales, abs(CWT_mor));  % Graficar la magnitud de la CWT
xlabel('Tiempo');
ylabel('Escalas');
title('CWT con ondita morlet (analitica)');
colorbar;
