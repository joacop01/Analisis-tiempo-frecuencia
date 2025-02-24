clear all;
close all;
addpath('../Funciones');

Fs = 1000;           % Frecuencia de muestreo (Hz)
N = 1000;            % Número de muestras
t = (0:N-1) / Fs;    % Vector de tiempo
f = 0 : Fs/N : Fs/2-Fs/N;  % Vector de frecuencias
D = 1;% Número de armónicos
f0 = 50;

n = 8;

x = exp(1i*2*pi*f0*t);
k = 200;

F = STFT_Gauss(x, t, k);

cant_crestas = 1;
Q = 5;

[c, ~] = Deteccion_Crestas(F, n, N, cant_crestas, Q);  

b = 20;

y = Reconstruccion_Cresta(c, b, N, F, cant_crestas);
y = y / (cant_crestas*N);

modulo = abs(y);
fase = angle(y)/(2*pi);
% Construir la matriz C
C = zeros(2*D, N);           % Inicialización de C

for d = 1:D
    C(d, :) = modulo.*cos(2 * pi *d * fase);    % Componentes coseno (c_d)
    C(D+d, :) = modulo.*sin(2 * pi * d * fase); % Componentes seno (d_d)
end

% Resolver el problema de optimización para obtener el vector v
v_estimado = (x * C') / (C * C');       % v = (xC^T)(CC^T)^(-1)
x_estimado = (v_estimado * C) + mean(x);


wave_func = Wave_function(v_estimado, D, length(x_estimado), Fs);

figure;
subplot(221);
plot(t, real(x_estimado));
hold on;
plot(t, real(x));
legend('x estimado', 'x');
title('Comparación de parte real');
hold off;
xlabel('Tiempo (s)');
ylabel('Amplitud');

subplot(223);
plot(t, imag(x_estimado));
hold on;
plot(t, imag(x));
legend('x estimado', 'x');
title('Comparación de parte imaginaria');
hold off;
xlabel('Tiempo (s)');
ylabel('Amplitud');

subplot(2, 2, [2, 4]);
plot(t, real(wave_func));
title('Función de onda');
xlabel('Tiempo (s)');
ylabel('Amplitud');
% 
F_estimado = STFT_Gauss(x_estimado, t, k);

figure;
subplot(211);
Plot_STFT(F_estimado, t, f);
title('STFT de la señal estimada');
subplot(212);
Plot_STFT(F, t, f);
title('STFT de la señal original');

%%
Fs = 1000;           % Frecuencia de muestreo (Hz)
N = 1000;            % Número de muestras
t = (0:N-1) / Fs;    % Vector de tiempo
f = 0 : Fs/N : Fs/2-Fs/N;  % Vector de frecuencias
D = 2;% Número de armónicos
f0 = 50;

n = 8;
indices = round(linspace(0.5/(n+1), 1 - 0.5/(n+1), n) * N);

x = exp(1i*2*pi*f0*t);
k = 200;

F = STFT_Gauss(x, t, k);

desvio_std =median(abs(real(F(:))))/0.6745;
Q_est = sqrt(2*gammaincinv(.997,1));
espectrograma = abs(F).^2;
umbral = Q_est * desvio_std;
F_somb = F;
F_somb (espectrograma < umbral) = 0;

cant_crestas = 1;
Q = 5;

[c, ~] = Deteccion_Crestas(F_somb, n, N, cant_crestas, Q);  

y = size(x);
b = 20;

for i =1:N

    I = c(i, 1) - b : c(i, 1) + b;

    I = max(I(1), 1):min(I(end), N/2);

    y(i) = sum(F_somb(I,i));
    y(i) = y(i) / N;
    
end

modulo = abs(y);
fase = angle(y)/(2*pi);
% Construir la matriz C

C = zeros(2*D, N);           % Inicialización de C

for d = 1:D
    C(d, :) = modulo.*cos(2 * pi *d * fase);    % Componentes coseno (c_d)
    C(d+D, :) = modulo.*sin(2 * pi * d * fase); % Componentes seno (d_d)
end

% Resolver el problema de optimización para obtener el vector v
% v = (y * C') / (C * C');       % v = (xC^T)(CC^T)^(-1)
media = mean(x);
x = x - mean(x);
[x_estimado, v_estimado] = Wave_Shape(y, D, N, x, mean(x));

wave_func = Wave_function(v_estimado, D, N, Fs);

figure;
subplot(221);
plot(t, real(x_estimado));
hold on;
plot(t, real(x));
legend('x estimado', 'x');
title('Comparación de parte real');
hold off;
xlabel('Tiempo (s)');
ylabel('Amplitud');

subplot(223);
plot(t, imag(x_estimado));
hold on;
plot(t, imag(x));
legend('x estimado', 'x');
title('Comparación de parte imaginaria');
hold off;
xlabel('Tiempo (s)');
ylabel('Amplitud');

subplot(2, 2, [2, 4]);
plot(t, real(wave_func));
title('Comparación de las formas de onda');
xlabel('Tiempo (s)');
ylabel('Amplitud');

F_estimado = STFT_Gauss(x_estimado, t, k);

figure;
subplot(211);
Plot_STFT(F_estimado, t, f);
title('STFT umbralada de la señal original');
subplot(212);
Plot_STFT(F_somb, t, f);
title('STFT de la señal estimada');

%%
Fs = 1000;           % Frecuencia de muestreo (Hz)
N = 1000;            % Número de muestras
t = (0:N-1) / Fs;    % Vector de tiempo
f = 0 : Fs/N : Fs/2-Fs/N;  % Vector de frecuencias
D = 2;
f0 = 50;

n = 8;

x = exp(1i*2*pi*f0*t);
k = 200;

F = STFT_Gauss(x, t, k);
F_g = STFT_Gauss_diff(x, t, k);
T = Synchro_STFT(F, F_g, t, N);

cant_crestas = 1;
Q = 5;

[c, ~] = Deteccion_Crestas(T, n, N, cant_crestas, Q);  

b = 20;

y = Reconstruccion_Cresta(c, b, N, F, cant_crestas);
y = y / N;
modulo = abs(y);
fase = angle(y)/(2*pi);

% Construir la matriz C
C = zeros(2*D, N);           % Inicialización de C
for d = 1:D
    C(d, :) = modulo.*cos(2 * pi *d * fase);    % Componentes coseno (c_d)
    C(d+D, :) = modulo.*sin(2 * pi * d * fase); % Componentes seno (d_d)
end

% Resolver el problema de optimización para obtener el vector v
% v = (y * C') / (C * C');       % v = (xC^T)(CC^T)^(-1)

media = mean(x);
x = x - media;
[x_estimado, v_estimado] = Wave_Shape(y, D, N, x, media);

wave_func = Wave_function(v_estimado, D, N, Fs);

figure;
subplot(221);
plot(t, real(x_estimado));
hold on;
plot(t, real(x));
legend('x estimado', 'x');
title('Comparación de parte real');
hold off;
xlabel('Tiempo (s)');
ylabel('Amplitud');

subplot(223);
plot(t, imag(x_estimado));
hold on;
plot(t, imag(x));
legend('x estimado', 'x');
title('Comparación de parte imaginaria');
hold off;
xlabel('Tiempo (s)');
ylabel('Amplitud');

subplot(2, 2, [2, 4]);
title('Comparación de las formas de onda');
plot(t, real(wave_func));
xlabel('Tiempo (s)');
ylabel('Amplitud');

F_estimado = STFT_Gauss(x_estimado, t, k);
figure;
subplot(212);
Plot_STFT(F_estimado, t, f);
title('STFT de la funcion estimada');
subplot(211);
Plot_STFT(T, t, f);
title('STFT con synchrosqueezing de la función original');
