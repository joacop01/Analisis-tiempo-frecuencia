clear all;
close all;
addpath('../Funciones');

f0 = 2;
f1 = 10;
Fs = 1000;
N = 1000;
Ts = 1 / Fs;
t = 0: Ts : N*Ts - Ts;
x1 = cos (2*pi*f0*t);
x2 = cos(2*pi*f1*t*2+1);
senial = x1 + x2;

modos = 2;

% 
% r = senial;
% IMFs = [];
% j = 1;
% k = 1;
% 
% while(k < modos)
%     while (true)  
%         [env_max, env_min] = encontrar_envolventes(r);
% 
%         m = (env_max + env_min) / 2;
% 
%         d = r - m;
% 
%         if(es_IMF(d))
%             IMFs = [IMFs, d];
%             r = r - d;
%             k = k +1;
%             j = 1;
%             break;
%         else
%             r = d;
%             j = j + 1;
%         end
%     end
%     
% end

% Subplot 1: Aproximación

[IMFs] = EMD(senial, modos);

[env_max, env_min] = encontrar_envolventes(senial);

figure;
plot(t, senial);
hold on;
plot(t, env_max);
plot(t, env_min);

figure;
subplot(321);
plot(t, IMFs(1,:));
title('IMF 1');
xlabel('Tiempo (s)');
ylabel('Amplitud');
xlim([0, 1]); % Limitar eje x entre 0 y 1 segundos
ylim([-2, 2]); % Limitar eje y entre -2 y 2

% Subplot 2: Detalle
subplot(322);
plot(t, IMFs(2, :));
title('IMF 2');
xlabel('Tiempo (s)');
ylabel('Amplitud');
xlim([0, 1]); % Limitar eje x entre 0 y 1 segundos
ylim([-2, 2]); % Limitar eje y entre -2 y 2

% Subplot 3: Señal de 50 Hz
subplot(323);
plot(t, x2);
title('Señal de 5 Hz');
xlabel('Tiempo (s)');
ylabel('Amplitud');
xlim([0, 1]); % Limitar eje x entre 0 y 1 segundos
ylim([-2, 2]); % Limitar eje y entre -2 y 2

% Subplot 4: Señal de 100 Hz
subplot(324);
plot(t, x1);
title('Señal de 2 Hz');
xlabel('Tiempo (s)');
ylabel('Amplitud');
xlim([0, 1]); % Limitar eje x entre 0 y 1 segundos
ylim([-2, 2]); % Limitar eje y entre -2 y 2

subplot(3, 2, [5 6]);
plot(t, senial);
title('Señal');
xlabel('Tiempo (s)');
ylabel('Amplitud');
xlim([0, 1]); % Limitar eje x entre 0 y 1 segundos
ylim([-2, 2]); % Limitar eje y entre -2 y 2


% Título general
sgtitle('Gráficos de Señales');















        