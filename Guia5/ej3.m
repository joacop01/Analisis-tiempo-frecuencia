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
f = -Fs/2: Fs/N : Fs/2 - Fs/N;

s = linspace(0,1,N);
s = transpose(s);
s_2 = repmat(s, 1, N);
% 
% x1 = cos(2*pi*100*t + 2*pi*100*t.^2);
% x2 = cos(2*pi*150*t + 2*pi*100*t.^2);
% x = x1 + x2;
% a = 20;
% b = 80;
% x = exp(1i*2*pi*(a.*t+b*t.^2));

x = exp(1i*2*pi*f0*t);

% Cálculo de la STFT y su derivada con ventana Gaussiana
CWT_mor = CWT_morlet(x, f, 180, linspace(0,1,N), 50);
CWT_mor_diff = CWT_morlet_diff(x, f, 180, linspace(0,1,N), 50);

delta_s = s(2) - s(1);


% Cálculo de la frecuencia instantánea
f_inst =   -imag((CWT_mor_diff./(CWT_mor*s_2*2*pi)));  % Frecuencia instantánea
% Inicializar la matriz de Synchrosqueezing
T = zeros(size(CWT_mor));

for n = 1: length(t)
    for k = 1:length(s)
        % Sincronización: reasignación de la energía
        a = f_inst(:, n);
        a = round(a);
        
        % Verificar si 'a' es un índice válido
        if a(k) >= 1 && a(k) <= length(s)
            if a(k) == s(k)
                T(a(k), n) = T(a(k), n) +  CWT_mor(k, n);
            end
        end
    end
end


%Graficar la CWT
figure;
subplot(311);
imagesc(t, linspace(0,N), abs(CWT_mor));  % Graficar la magnitud de la CWT
xlabel('Tiempo');
ylabel('Escalas');
title('CWT con ondita morlet (analitica)');
colorbar;

subplot(312);
imagesc(t, linspace(0,N), abs(CWT_mor_diff));  % Graficar la magnitud de la CWT
xlabel('Tiempo');
ylabel('Escalas');
title('CWT con ondita morlet derivada (analitica)');
colorbar;

% Graficar la representación Synchrosqueezed
subplot(313);
imagesc(t, linspace(0,N), abs(T));  % Graficar la magnitud de la CWT
xlabel('Tiempo');
ylabel('Escalas');
title('Synchrosqueezed CWT');
colorbar;

