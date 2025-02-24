clear all;
close all;
addpath('../Funciones');
%Mexhat
% Parámetros de muestreo
Fs = 1000;             % Frecuencia de muestreo
Ts = 1/Fs;             % Periodo de muestreo
t = 0:Ts:1-Ts;             % Vector de tiempo
N = length(t);
f = -Fs/2:Fs/N: Fs/2-Fs/N; %Vector de frecuencias 
k = 0:Fs/N: Fs/2-Fs/N; 

% Señal de prueba
% x = cos(2*pi*100*t);     % Ejemplo de señal senoidal
x1 = cos(2*pi*100*t + 2*pi*100*t.^2);
x2 = cos(2*pi*150*t + 2*pi*100*t.^2);

x = x1 ;%+ x2;
figure;
plot(t, x);


% Parámetros de la wavelet sombrero mexicano
sigma = 500;           % Parámetro de escala de la Gaussiana

F = STFT_Gauss(x, t, sigma);
figure;
Plot_STFT(F, t, k);

% FFT de la señal
x_f = fftshift(fft(x));

% Definir las escalas (puedes ajustar los valores de escala)
scales = 0.001:0.0001:0.5;
% Inicializar la matriz de CWT
CWT = zeros(length(scales),length(x));
    
for k = 1: length(scales)
    
    f_scaled = (f*scales(k));
    mex_f = 2*pi^2*(f_scaled.^2)*sqrt(pi/sigma).*exp(-pi^2*(f_scaled.^2)/(sigma));


    % Multiplicación en frecuencia (señal * conj(wavelet escalada))
    X_wav_f = x_f .* conj(mex_f);

    % Inversa de Fourier para obtener la CWT en el dominio del tiempo
    CWT(k, :) = ifft(X_wav_f);
end

% Inicializar la matriz de CWT
CWT_analitica = CWT_mexhat(x, f, sigma, scales);

%Graficar la CWT
figure;
subplot(211);
imagesc(t, scales,abs(CWT_analitica));  % Graficar la magnitud de la CWT
xlabel('Tiempo');
ylabel('Escalas');
title('CWT con ondita mexhat (analitica)');
colorbar;
subplot(212);
imagesc(t, scales, abs(CWT));  % Graficar la magnitud de la CWT
xlabel('Tiempo');
ylabel('Escalas');
title('CWT con ondita mexhat (real)');
colorbar;
