clear all;
close all;
addpath('../Funciones');

Fs = 1000;           % Frecuencia de muestreo (Hz)
N = 1000;            % Número de muestras

t = (0:N-1) / Fs;    % Vector de tiempo
f = 0 : Fs/N : Fs/2-Fs/N;  % Vector de frecuencias
D = 1;% Número de armónicos
a = 50;
b = 40;
modulo = 10;

[x_generado, v] = Sintetic_cuad_func(a, b, D, modulo, t);

figure;
subplot(311);
plot(t, x_generado);
title('Señal generada');
xlabel('Tiempo (s)');
ylabel('Amplitud');

F = STFT_Gauss(x_generado, t, 500);

subplot(312);
Plot_STFT(F, t, f);

wave_func = Wave_function(v, D, N, Fs);

subplot(313);
plot(t, wave_func);
title('Forma de onda');
sgtitle('Fase cuadrática');
xlabel('Tiempo (s)');
ylabel('Amplitud');

%%
Fs = 1000;           % Frecuencia de muestreo (Hz)
N = 1000;            % Número de muestras
t = (0:N-1) / Fs;    % Vector de tiempo
f = 0 : Fs/N : Fs/2-Fs/N;  % Vector de frecuencias
D = 1;% Número de armónicos
modulo = 10;

[x_generado, v] = Sintetic_cos_func(a, b, D, modulo, t);

figure;
subplot(311);
plot(t, x_generado);
title('Señal generada');
xlabel('Tiempo (s)');
ylabel('Amplitud');

F = STFT_Gauss(x_generado, t, 700);

subplot(312);
Plot_STFT(F, t, f);

wave_func = Wave_function(v, D, N, Fs);

subplot(313);
plot(t, wave_func);
title('Forma de onda');
sgtitle('Fase cosenoidal');
xlabel('Tiempo (s)');
ylabel('Amplitud');