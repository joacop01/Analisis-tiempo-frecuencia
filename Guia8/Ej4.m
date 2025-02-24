clear all;
close all;
addpath('../Funciones');



x1 = load('rr1.txt');
x1 = transpose(x1);
x2 = load('rr2.txt');
x2 = transpose(x2);

N1 = length(x1);
N2 = length(x2);

t1 = 1 : N1;
t2 = 1 : N2;

figure;
subplot(211);
plot(t1, x1);
title('rr1');

subplot(212);
plot(t2, x2);
title('rr2')

modos = 8;

IMFs1 = EMD(x1, modos);
IMFs2 = EMD(x2, modos);

figure;
for j = 1: modos
    subplot(3, 4, j);
    plot(t2, IMFs2(j, :));
    title(sprintf('IMF %d', j)); % Cambia dinámicamente el número de IMF
    xlabel('Tiempo (s)');
    ylabel('Amplitud');

end
subplot(3, 4, [9 10 11 12]);
plot(t2, x2);
title('Señal');
xlabel('Tiempo (s)');
ylabel('Amplitud');
sgtitle('Descomposición empírica en modos');

tendencia2 = zeros(1, N2);

for k = 2: modos
    tendencia2 = tendencia2 + IMFs2(k, :);
end

tendencia1 = zeros(1, N1);

for k = 2: modos
    tendencia1 = tendencia1 + IMFs1(k, :);
end

figure;
subplot(211);
plot(t1, x1);
hold on;
plot(t1, tendencia1);
title('Señal y su tendencia');
legend('rr1', 'tendencia');
hold off;

subplot(212);
plot(t2, x2);
hold on;
plot(t2, tendencia2);
title('Señal y su tendencia');
legend('rr2', 'tendencia');
hold off;