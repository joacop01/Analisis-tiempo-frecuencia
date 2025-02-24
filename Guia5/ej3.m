clear all;
close all;
addpath('../Funciones');
%Parámetros de la señal
tol = 1;
Fs = 1000;       % Frecuencia de muestreo
Ts = 1/Fs;       % Periodo de muestreo
t = 0:Ts:1-Ts;   % Eje temporal
k = 500;         % Parámetro de ventana
f0 = 200;        % Frecuencia de la señal
N = length(t);   % Número de puntos
f = -Fs/2: Fs/N : Fs/2 - Fs/N;

s = linspace(0.01,1,N);
s = transpose(s);
s_2 = repmat(s, 1, N);

x = exp(1i*2*pi*f0*t);

% Cálculo de la STFT y su derivada con ventana Gaussiana
CWT_mor = CWT_morlet(x,f, 200, s, 100);
CWT_mor_diff = CWT_morlet_diff(x,f, 200, s, 100);

delta_s = s(2) - s(1);

% Cálculo de la frecuencia instantánea
f_inst =   real(-CWT_mor_diff./(CWT_mor.*s_2*2*pi*1j));  % Frecuencia instantánea

% Inicializar la matriz de Synchrosqueezing
T = zeros(size(CWT_mor));

for n = 1: length(t)
    for k = 1:length(s)
        % Sincronización: reasignación de la energía
        a = f_inst(:, n);
%         Verificar si 'a' es un índice válido
        if abs(a(k) - s(k)) < tol
        
            T(k, n) = T(k, n) + CWT_mor(k, n);
        end    
    end
end


%Graficar la CWT
figure;
subplot(311);
imagesc(t, s, abs(CWT_mor));  % Graficar la magnitud de la CWT
xlabel('Tiempo');
ylabel('Escalas');
title('CWT con ondita morlet (analitica)');
colorbar;

subplot(312);
imagesc(t, s, abs(CWT_mor_diff));  % Graficar la magnitud de la CWT
xlabel('Tiempo');
ylabel('Escalas');
title('CWT con ondita morlet derivada (analitica)');


% Graficar la representación Synchrosqueezed
subplot(313);
Plot_STFT(T, t, f);
title('Synchrosqueezed CWT');

