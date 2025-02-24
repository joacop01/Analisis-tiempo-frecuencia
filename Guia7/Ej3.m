clear all;
close all;
addpath('../Funciones');

Fs = 1000;           % Frecuencia de muestreo (Hz)
N = 1000;            % Número de muestras
Ts = 1/Fs;
t = (0:N-1) / Fs;
f = 0 : Fs/N : Fs/2-Fs/N;  % Vector de frecuencias

a = 150;
b = 100;
D = 2;
modulo = 10;

[cos_sint, v_cos] = Sintetic_cos_func(a, b, D, modulo, t);

b = 40;

[cuad_sint, v_cuad] = Sintetic_cuad_func(a, b, D, modulo, t);

snr_db = 10;
P_signal = mean(cos_sint.^2);

P_noise = P_signal / 10^(snr_db / 10);

% Generar ruido gaussiano con la potencia deseada
ruido_gaussiano = sqrt(P_noise) * randn(size(cos_sint));  % Escalar el ruido

% Añadir ruido a la señal y almacenar en la columna correspondiente
x_ruido = cos_sint + ruido_gaussiano;  % Transponer para que cada columna corresponda a una señal ruidosa
P_noise = mean(ruido_gaussiano.^2);  % Potencia del ruido
SNR = 10 * log10(P_signal / P_noise);

% Mostrar el valor de la SNR (debería ser cercano a 0 dB)
disp(['Relación señal-ruido (SNR): ', num2str(SNR), ' dB']);


figure;
subplot(411);
plot(t, cos_sint);
title('Señal original');
xlabel('Tiempo (s)');
ylabel('Amplitud');

subplot(412);
plot(t, x_ruido);
title('Señal contaminada');
xlabel('Tiempo (s)');
ylabel('Amplitud');

k = 1000;
F = STFT_Gauss(cos_sint, t, k);

F_ruido = STFT_Gauss(x_ruido, t, k);


subplot(413);
Plot_STFT(F, t, f);
subplot(414);
Plot_STFT(F_ruido, t, f);
title('STFT de la señal ruidosa');

n = 8;
cant_crestas = 1;
Q = 5;

[c, ~] = Deteccion_Crestas(F_ruido, n, N, cant_crestas, Q);  

y = Reconstruccion_Cresta(c, b, N, F, cant_crestas);
y = y / (N/2);

media = mean(x_ruido);
x_ruido = x_ruido - media;
[x_estimado, v_estimado] = Wave_Shape(y, D, N, x_ruido ,media);

[wave_func, t_wave] = Wave_function(v_estimado, D, N, Fs);
[wave_func_orig, t_wave_orig] = Wave_function(v_cos, D, N, Fs);

figure;
subplot(2, 2, [1 2]);
plot(t, cos_sint);
hold on;
plot(t, x_estimado);
title('Función estimada');
xlabel('Tiempo (s)');
ylabel('Amplitud');
hold off;

subplot(2, 2, [3 4]);
plot(t_wave_orig, wave_func_orig);
hold on;
plot(t_wave, wave_func);
title('Funciones de onda');
xlabel('Tiempo (s)');
ylabel('Amplitud');
legend('Función de onda original', 'Función de onda estimada');
hold off;

sgtitle('Comparación de las funciones');

F_estimado = STFT_Gauss(x_estimado, t, k);

figure;
subplot(211);
Plot_STFT(F_ruido, t, f);
subplot(212);
Plot_STFT(F_estimado, t, f);
title('STFT de la reconstrucción');
