addpath('C:\Repositorios\Curso Análisis Tiempo-Frecuencia\Funciones');
%Parámetros de la señal
Fs = 1000;       % Frecuencia de muestreo
Ts = 1/Fs;       % Periodo de muestreo
t = 0:Ts:1-Ts;   % Eje temporal
k = 200;         % Parámetro de ventana
f0 = 300;        % Frecuencia de la señal
N = length(t);   % Número de puntos

x1 = cos(2*pi*100*t + 2*pi*100*t.^2);
x2 = cos(2*pi*150*t + 2*pi*100*t.^2);
x = x1 + x2;

% Cálculo de la STFT y su derivada con ventana Gaussiana
F = STFT_Gauss(x, t, k);        % STFT
F_g = STFT_Gauss_diff(x, t, k); % Derivada de la STFT

% Eje de frecuencias
f = 0: Fs/N : Fs/2 - Fs/N;

% Cálculo de la frecuencia instantánea
f_inst = (1/(2*pi))*imag(F_g./ F);  % Frecuencia instantánea
% Inicializar la matriz de Synchrosqueezing
T = zeros(size(F));
a = zeros(size(F));
for n = 1:length(t)
    for m = 1:length(f)
        % Sincronización: reasignación de la energía
        d = 100000*abs(f_inst(m, n));
        d = round(d);
        a(m,n) = d; 
        
        % Verificar si 'a' es un índice válido
        if a(m,n) >= 1 && a(m,n) <= length(f)
            if a == f(m)
                T(a(m,n), n) = T(a(m,n), n) + sum(F(m, :));
            end
        end
    end
end

% Bucle sobre el tiempo y la frecuencia
% Graficar la representación Synchrosqueezed
figure;
imagesc(t, f, abs(F(Fs/2:end, :)));
axis xy;
xlabel('Tiempo (s)');
ylabel('Frecuencia (Hz)');
title('Synchrosqueezed STFT');
colorbar;

% Graficar la representación Synchrosqueezed
figure;
imagesc(t, f, abs(T(Fs/2:end, :)));
axis xy;
xlabel('Tiempo (s)');
ylabel('Frecuencia (Hz)');
title('Synchrosqueezed STFT');
colorbar;