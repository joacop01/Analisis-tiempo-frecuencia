clear all;
close all;
addpath('../Funciones');

Fs = 1000;             % Frecuencia de muestreo
Ts = 1/Fs;             % Periodo de muestreo
t = 0:Ts:1-Ts;         % Vector de tiempo
N = length(t);         % Número de puntos en el tiempo

f = 0 : Fs/N : Fs/2-Fs/N;  % Vector de frecuencias
Q = 5;

x1 = cos(2*pi*100*t + 2*pi*100*t.^2);
x2 = cos(2*pi*150*t + 2*pi*100*t.^2);

x = x1 + x2;
% x = cos(2*pi*100*t);

F = STFT_Gauss(x, t, 200);
n = 8;
indices = round(linspace(0.5/(n+1), 1 - 0.5/(n+1), n) * N);

c = zeros(N,1);

[~, c(indices(1))] = max(abs(F(:,indices(1)).^2));

for a = indices(1)+1:N
    I_a = c(a-1) - Q: c(a-1) + Q;
    [~, indice_local_a] = max(abs(F(I_a, a)).^2);
    
    c(a) = I_a(indice_local_a);

end

for b = indices(1)-1:-1:1
    I_b = c(b+1) - Q: c(b+1) + Q;
    [~, indice_local_b] = max(abs(F(I_b, b)).^2);
    
    c(b) = I_b(indice_local_b);
end

%Borramos la cresta recien detectada
F_aux = F;
for k = 1:N
    I = c(k) - Q : c(k) + Q;
    F_aux(I, k) = 0;
end

%Deteccion de la segunda cresta
c_2 = zeros(N,1);

[~, c_2(indices(1))] = max(abs(F_aux(:,indices(1)).^2));
for a = indices(1)+1:N
    I_a = c_2(a-1) - Q: c_2(a-1) + Q;
    [~, indice_local_a] = max(abs(F_aux(I_a, a)).^2);
    
    c_2(a) = I_a(indice_local_a);

end

for b = indices(1)-1:-1:1
    I_b = c_2(b+1) - Q: c_2(b+1) + Q;
    [~, indice_local_b] = max(abs(F_aux(I_b, b)).^2);
    
    c_2(b) = I_b(indice_local_b);
end



Plot_STFT(F, t, f);
hold on;
plot(t, c, 'r');
plot(t, c_2, 'b');
legend('Cresta 1','Cresta 2');
hold off;
