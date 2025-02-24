clear all;
close all;
addpath('../Funciones');

N = 1000;
n = 1: N;

x1 = sin(2*pi*(10^-5)*n.^2);

x2 = sin(2*pi*(10^-5)*n.^2 + 2*pi*(10^-2)*n);

x = x1 + x2;

figure;
subplot(311);
plot(n, x);

subplot(312);
plot(n, x1);

subplot(313);
plot(n, x2);

modos = 2;
IMFs = EMD(x, modos);

[env_max, env_min] = encontrar_envolventes(x);

figure;
plot(n, x);
hold on;
plot(n, env_max);
plot(n, env_min);
plot(n, (env_max + env_min)/2);
hold off;


figure;
subplot(321);
plot(n,IMFs(1, :));
title('IMF 1');
xlabel('Tiempo (s)');
ylabel('Amplitud');

subplot(322);
plot(n, IMFs(2, :));
title('IMF 2');
xlabel('Tiempo (s)');
ylabel('Amplitud');

subplot(323);
plot(n, x2);
title('Modo 1');
xlabel('Tiempo (s)');
ylabel('Amplitud');

subplot(324);
plot(n, x1);
title('Modo 2');
xlabel('Tiempo (s)');
ylabel('Amplitud');

subplot(3, 2, [5 6]);
plot(n, x);
title('Señal');
xlabel('Tiempo (s)');
ylabel('Amplitud');

% Título general
sgtitle('Gráficos de Señales');


% IMFs = zeros(2, N);
% k = 1;
% r = x;
% cantidad = 7;
% for j = 1 : cantidad  
%     [env_max, env_min] = encontrar_envolventes(r);
% 
%     m = (env_max + env_min) / 2;
% 
%     d = r - m;
%     if (j == cantidad-1 || j == 1)
%         figure;
%         subplot(211);
%         plot(n, r);
%         hold on;
%         plot(n, env_max);
%         plot(n, env_min);
%         plot(n, m);
%         subplot(212);
%         plot(n, d);
%     end
%     if(es_IMF(d))
%         IMFs(k, :) = d;
%         r = m;
%         k = k +1;
%         disp(num2str(j));
%         disp('La señal es IMF');
%         figure;
%         subplot(211);
%         plot(n, d);
%         subplot(212);
%         plot(n, r);
%         break;
%     else
%     r = d;
% %         j = j + 1;
%     end
% end