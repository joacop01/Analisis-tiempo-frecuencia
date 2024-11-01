close all;
addpath('../Funciones');
%Parámetros de la señal

Fs = 1000;       % Frecuencia de muestreo
Ts = 1/Fs;       % Periodo de muestreo
t = 0:Ts:1-Ts;   % Eje temporal
k = 1500;         % Parámetro de ventana
f0 = 1;        % Frecuencia de la señal
N = length(t);   % Número de puntos
f = 0: Fs/N : Fs/2 - Fs/N;

% x1 = cos(2*pi*100*t + 2*pi*100*t.^2);
% x2 = cos(2*pi*150*t + 2*pi*100*t.^2);
% x = x1 + x2;
x1 = cos(2*pi*(150*t+ 100/(2*pi)*sin(2*pi*t)));
x2 = cos(2*pi*(300*t+120/(2*pi)*sin(2*pi*t)));

x = x1 + x2;

% a = 100;
% phi = sin(2*pi*f0*t); 
% b = 10;
% x = exp(1i*2*pi*(a.*t+b*phi));

% Calcular la potencia de la señal
P_signal = mean(x.^2);  % Potencia de la señal original

% Generar ruido gaussiano con la misma potencia que la señal
ruido_gaussiano = sqrt(P_signal) * randn(size(x));  % Escalar el ruido para que tenga la misma potencia

%x = x + ruido_gaussiano;


% x = exp(1i*2*pi*f0*t);

F = STFT_Gauss(x, t, k);
F_g = STFT_Gauss_diff(x, t, k);
T = Synchro_STFT(F, F_g, t, N);

% Graficar la STFT
Plot_STFT(F, t, f);
title("|STFT|");
Plot_STFT(T, t, f);
title('|Synchorsqueezed STFT|');

% Graficar la representación Synchrosqueezed