clear all;
close all;
addpath('../Funciones');

Fs = 100;
Ts = 1/Fs;
N = 1000;
t = 1: N;
k = 10;

x1 = sin(2*pi*0.1*t);
v_g = exp(-k*((t-500)/1000).^2);
x1 = x1.*v_g;

% figure;
% subplot(211);
% plot(t, x1);
% subplot(212);
% plot(t, v_g);

t_mod = mod(t, 400);
x2 = (t_mod >= 0 & t_mod <= 200) .* (t_mod / 100) + (t_mod > 200 & t_mod <= 400) .* (4 - (t_mod/ 100));


x = x1 + x2;

figure;
subplot(311);
plot(t, x);
subplot(312);
plot(t, x1);
subplot(313);
plot(t, x2);


modos = 2;
[env_max, env_min] = encontrar_envolventes(x);

figure;
plot(t, x);
hold on;
plot(t, env_max);
plot(t, env_min);

IMFs = EMD(x, modos);

figure;
subplot(221);
plot(t,IMFs(1, :));
title('IMF 1');
xlabel('Tiempo (s)');
ylabel('Amplitud');
ylim([-1, 1]); % Limitar eje y entre -2 y 2

subplot(222);
plot(t, IMFs(2, :));
title('IMF 2');
xlabel('Tiempo (s)');
ylabel('Amplitud');

subplot(223);
plot(t, x1);
title('Modo 1');
xlabel('Tiempo (s)');
ylabel('Amplitud');

subplot(224);
plot(t, x2);
title('Modo 2');
xlabel('Tiempo (s)');
ylabel('Amplitud');

% Título general
sgtitle('Gráficos de Señales');



% r = x;
% cantidad = 1;
% for j = 1 : cantidad  
%     [env_max, env_min] = encontrar_envolventes(r);
% 
%     m = (env_max + env_min) / 2;
% 
%     d = r - m;
%     if (j == cantidad-1 || j == 1)
%         figure;
%         subplot(211);
%         plot(t, r);
%         hold on;
%         plot(t, env_max);
%         plot(t, env_min);
%         plot(t, m);
%         subplot(212);
%         plot(t, d);
%     end
%     if(es_IMF(d))
%         IMFs(k, :) = d;
%         r = r- d;
%         k = k +1;
%         disp('La señal es IMF');
%         figure;
%         subplot(211);
%         plot(t, d);
%         subplot(212);
%         plot(t, r);
%         break;
%     else
%         r = d;
% %         j = j + 1;
%     end
% end





% 
% [pks_max, locs_max] = findpeaks(x, t); % Máximos locales
% [pks_min, locs_min] = findpeaks(-x, t); % Mínimos locales (invertir la señal)
% pks_min = -pks_min; % Revertir el signo de los mínimos
%  % Agregar los extremos iniciales y finales si no están incluidos
% if t(1) < locs_max(1) || x(1) > x(2)
%     locs_max = [t(1), locs_max];
%     pks_max = [x(1), pks_max];
% end
% if t(end) > locs_max(end) || x(end) > x(end-1)
%     locs_max = [locs_max, t(end)];
%     pks_max = [pks_max, x(end)];
% end
% 
% if t(1) < locs_min(1) || x(1) < x(2)
%     locs_min = [t(1), locs_min];
%     pks_min = [x(1), pks_min];
% end
% if t(end) > locs_min(end) || x(end) < x(end-1)
%     locs_min = [locs_min, t(end)];
%     pks_min = [pks_min, x(end)];
% end
% 
% env_max = interp1(locs_max, pks_max, t, 'pchip'); % Envolvente máxima
% env_min = interp1(locs_min, pks_min, t, 'pchip'); % Envolvente mínima
% 
%     % Graficar la señal, extremos y envolventes
% figure;
% plot(t, x, 'b', 'LineWidth', 1.2); 
% hold on;
% plot(locs_max, pks_max, 'ro', 'MarkerSize', 6, 'DisplayName', 'Máximos Ajustados');
% plot(locs_min, pks_min, 'go', 'MarkerSize', 6, 'DisplayName', 'Mínimos Ajustados');
% plot(t, env_max, 'r--', 'LineWidth', 1.2, 'DisplayName', 'Envolvente Máxima');
% plot(t, env_min, 'g--', 'LineWidth', 1.2, 'DisplayName', 'Envolvente Mínima');
% legend('Señal', 'Máximos Ajustados', 'Mínimos Ajustados', 'Envolvente Máxima', 'Envolvente Mínima');
% xlabel('Tiempo');
% ylabel('Amplitud');
% title('Detección Mejorada de Extremos y Envolventes');
% grid on;
% hold off;
