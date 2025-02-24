%%
%Morlet
clear all;
close all;
addpath('../Funciones');
Fs = 1000;             % Frecuencia de muestreo
Ts = 1/Fs;             % Periodo de muestreo
t = 0:Ts:1-Ts;         % Vector de tiempo
N = length(t);         % Número de puntos en el tiempo
f = -Fs/2 : Fs/N : Fs/2-Fs/N;  % Vector de frecuencias
sigma = 50;
f0 = 50;


x1 = cos(2*pi*100*t + 2*pi*100*t.^2);
x2 = cos(2*pi*150*t + 2*pi*100*t.^2);
x = x1 + x2;


% FFT de la señal
x_f = fftshift(fft(x));

% Definir las escalas (puedes ajustar los valores de escala)
scales = linspace(0,1,N);

% Inicializar la matriz de CWT
CWT = zeros(length(scales),length(x));

% Realizar la CWT en el dominio de la frecuencia
for k = 1:length(scales)
    
    morlet_f = sqrt(pi/sigma)*exp((-pi^2*(f*scales(k)-f0).^2)/sigma);
    
    % Multiplicación en frecuencia (señal * conj(wavelet escalada))
    X_wav_f = x_f .* conj(morlet_f);
    
    % Inversa de Fourier para obtener la CWT en el dominio del tiempo
    CWT(k, :) = ifft(X_wav_f);
end

CWT_analitica = CWT_morlet(x, f, 180, linspace(0,1,N), 50);

%Graficar la CWT
figure;
subplot(211);
imagesc(t, linspace(0,N), abs(CWT_analitica));  % Graficar la magnitud de la CWT
xlabel('Tiempo');
ylabel('Escalas');
title('CWT con ondita morlet(analitica)');
colorbar;
subplot(212);
imagesc(t, linspace(0,N), abs(CWT));  % Graficar la magnitud de la CWT
xlabel('Tiempo');
ylabel('Escalas');
title('CWT con ondita morlet (real)');
colorbar;

