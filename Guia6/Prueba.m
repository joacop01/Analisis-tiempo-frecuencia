Fs = 1000;             % Frecuencia de muestreo
Ts = 1/Fs;             % Periodo de muestreo
t = 0:Ts:1-Ts;         % Vector de tiempo
N = length(t);         % Número de puntos en el tiempo
cant_crestas = 1;
Q = 5;
n = 8;
indices = round(linspace(0.5/(n+1), 1 - 0.5/(n+1), n) * N);
f = 0 : Fs/N : Fs/2-Fs/N;  % Vector de frecuencias

% x = cos(2*pi*(100*t + 100*t.^2));
x = exp(i*2*pi*100*t);

F = STFT_Gauss(x, t, 100);
c = Deteccion_Crestas(F, indices, N, cant_crestas, Q);
Plot_STFT(F, t, f);
title(['Detección de cresta']);
hold on;
plot(t, c(:,1, 1), 'r');
plot(t, c, 'b');
legend('Cresta 1');
hold off;

y = size(x);
b = 18;
for i =1:N
    y(i) = sum((F(:,i)));
end

y = y./N;

figure;
plot(t, real(y));
hold on;
plot(t, imag(y));
plot(t, imag(x));
legend('Real', 'Imaginaria', 'Imaginaria original');
hold off;

% figure;
% % plot(t, imag(x));
% hold on;
% plot(t, imag(y));
% legend('Real', 'Reconstruccion');
% hold off;
