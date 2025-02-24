%%
%Bump
clear all;
close all;
addpath('../Funciones');
% Parámetros de muestreo
Fs = 1000;             % Frecuencia de muestreo
Ts = 1/Fs;             % Periodo de muestreo
t = 0:Ts:1-Ts;         % Vector de tiempo
N = length(t);         % Número de puntos en el tiempo
f = -Fs/ 2 : Fs/N : Fs/2-Fs/N;  % Vector de frecuencias
sigma = 8;             % Ancho de la protuberancia
u = 50;                 % Valor central

% Señal de prueba
x = cos(2*pi*200*t);     % Ejemplo de señal senoidal
% x1 = cos(2*pi*100*t + 2*pi*100*t.^2);
% x2 = cos(2*pi*150*t + 2*pi*100*t.^2);

% x = x1 + x2;
% FFT de la señal
x_f = fft(x);

% Definir las escalas (puedes ajustar los valores de escala)
scales = linspace(0,1,N);

% Inicializar la matriz de CWT
CWT = zeros(length(scales),length(x));

f_abs = abs(f);
% Realizar la CWT en el dominio de la frecuencia
for k = 1:length(scales)
        
    denominator = 1 - ((f*scales(k) - u).^2)/(sigma^2);  % Calcula el denominador de la función bump
    denominator(denominator == 0) = eps; % Asegura que no sea cero o negativo
    bump = exp(1 - (1 ./ (denominator))); % Función bump

    % Crear una envolvente en función del tiempo
    chi = (f*scales(k) >= u - sigma) & (f*scales(k) <= u + sigma); % Define la ventana temporal

    % Aplicar la envolvente
    bump = bump .* chi; % Aplicar la envolvente a la función bump

    X_wav_f = x_f .* conj(bump);

    % Inversa de Fourier para obtener la CWT en el dominio del tiempo
    CWT(k, :) = ifft(X_wav_f);
    
end

CWT_analitica = CWT_bump(x, f, 8, scales, 50);

%Graficar la CWT
figure;
subplot(211);
imagesc(t, linspace(0, N), abs(CWT_analitica));  % Graficar la magnitud de la CWT
xlabel('Tiempo');
ylabel('Escalas');
title('CWT con ondita bump (analitica)');
colorbar;
subplot(212);
imagesc(t, linspace(0, N), abs(CWT));  % Graficar la magnitud de la CWT
xlabel('Tiempo');
ylabel('Escalas');
title('CWT con ondita bump (real)');
colorbar;

