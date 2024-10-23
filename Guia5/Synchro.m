clear all;
close all;
addpath('../Funciones');
%Parámetros de la señal

Fs = 1000;       % Frecuencia de muestreo
Ts = 1/Fs;       % Periodo de muestreo
t = 0:Ts:1-Ts;   % Eje temporal
k = 300;         % Parámetro de ventana
f0 = 100;        % Frecuencia de la señal
N = length(t);   % Número de puntos
f = 0: Fs/N : Fs/2 - Fs/N;
f = transpose(f);
f_2 = repmat(f, 1, N);
% 
% x1 = cos(2*pi*100*t + 2*pi*100*t.^2);
% x2 = cos(2*pi*150*t + 2*pi*100*t.^2);
% x = x1 + x2;
a = 20;
b = 80;
x = exp(1i*2*pi*(a.*t+b*t.^2));

%x = exp(1i*2*pi*f0*t);

% Cálculo de la STFT y su derivada con ventana Gaussiana
F = STFT_Gauss(x, t, k);        % STFT
F_g = STFT_Gauss_diff(x, t, k); % Derivada de la STFT

% Cálculo de la frecuencia instantánea
f_inst =   f_2-(1/(2*pi)).*imag(F_g./ F);  % Frecuencia instantánea
% Inicializar la matriz de Synchrosqueezing
T = zeros(size(F));

for n = 1: length(t)
    for k = 1:length(f)
        % Sincronización: reasignación de la energía
        a = f_inst(:, n);
        a = round(a);
        
        % Verificar si 'a' es un índice válido

        if a(k) >= 1 && a(k) <= length(f)
            if a(k) == f(k)
                T(a(k), n) = T(a(k), n) +  (F(k, n));
            end
        end
        
    end  
end

figure;
subplot(211);
imagesc(t, f, abs(F));
axis xy;
xlabel('Tiempo (s)');
ylabel('Frecuencia (Hz)');
title('STFT');
colorbar;


% Graficar la representación Synchrosqueezed
subplot(212);
imagesc(t, f, abs(T));
axis xy;
xlabel('Tiempo (s)');
ylabel('Frecuencia (Hz)');
title('Synchrosqueezed STFT');
colorbar;