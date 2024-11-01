clear all;
close all;
addpath('../Funciones');

Fs = 1000;             % Frecuencia de muestreo
Ts = 1/Fs;             % Periodo de muestreo
t = 0:Ts:1-Ts;         % Vector de tiempo
N = length(t);         % Número de puntos en el tiempo
n = 8;
indices = round(linspace(0.5/(n+1), 1 - 0.5/(n+1), n) * N);
Q = 30;
f = 0 : Fs/N : Fs/2-Fs/N;  % Vector de frecuencias
x1 = cos(2*pi*(150*t+ (100/(2*pi))*sin(2*pi*t)));
x2 = cos(2*pi*(300*t+(120/(2*pi))*sin(2*pi*t)));

x = x1 + 0.25*x2;


% F = STFT_Gauss(x, t, 1800);
cant_crestas = 2;
% c = Deteccion_Crestas(F, indices, N, cant_crestas, Q);

y = real(exp(1i*2*pi*(150*t + 100/(2*pi)*sin(2*pi*t))) + 0.25*exp(1i*2*pi*(300*t + 120/(2*pi)*sin(2*pi*t))));
frec_inst_1 = 150 + 100*cos(2*pi*t);
frec_inst_2 = 300 + 120*cos(2*pi*t);

frec_inst_1 = transpose(frec_inst_1);
frec_inst_2 = transpose(frec_inst_2);

F_2 = STFT_Gauss(y, t, 1500);

Plot_STFT(F_2, t, f);
hold on;
plot(t, frec_inst_1,'r');
plot(t, frec_inst_2,'b');
title('Frecuencia instantánea real');
legend('IF1', 'IF2');
hold off;


P_signal = mean(x.^2);                % Potencia de la señal original
% Escalar el ruido para diferentes SNR
% Escalar el ruido para diferentes SNR
SNR_values = [0, 10, 20, 30];            % Valores de SNR en dB
ruido_gaussiano = zeros(size(x));     % Inicializar el ruido
x_ruido = zeros(length(SNR_values), length(x));  % Inicializar matriz para las señales ruidosas
P_noise = zeros(length(SNR_values), 1);

crestas_1 = zeros(N, 10);
crestas_2 = zeros(N, 10);

% for w = 1: 10
    
    for idx = 1:length(SNR_values)
        snr_db = SNR_values(idx);          % Obtener el valor de SNR actual

        % Calcular la potencia del ruido deseada
        P_noise(idx) = P_signal / 10^(snr_db / 10);

        % Generar ruido gaussiano con la potencia deseada
        ruido_gaussiano = sqrt(P_noise(idx)) * randn(size(x));  % Escalar el ruido

        % Añadir ruido a la señal y almacenar en la columna correspondiente
        x_ruido(idx, :) = x + ruido_gaussiano;  % Transponer para que cada columna corresponda a una señal ruidosa
        P_noise(idx) = mean(ruido_gaussiano.^2);  % Potencia del ruido
        SNR = 10 * log10(P_signal / P_noise(idx));

        % Mostrar el valor de la SNR (debería ser cercano a 0 dB)
        disp(['Relación señal-ruido (SNR): ', num2str(SNR), ' dB']);
    end

    for u = 1:length(SNR_values)
        F = STFT_Gauss(x_ruido(u,:), t, 1500);
        c = Deteccion_Crestas(F, indices, N, cant_crestas, Q);
        Plot_STFT(F, t, f);
        title(['Detección de cresta con ' num2str(SNR_values(u)) 'dB de ruido']);
        hold on;
        plot(t, c(:,1, 1), 'r');
        plot(t, c(:,2, 1), 'b');
        legend('Cresta 1','Cresta 2');
        hold off;


        ECT1 = sum((abs(c(:,1, 1)-frec_inst_1)).^2)/N;
        ECT2 = sum((abs(c(:,2, 1)-frec_inst_2)).^2)/N;
        disp(['ECT1: ', num2str(ECT1), ' para ', num2str(SNR_values(u)), ' dB']);
        disp(['ECT2: ', num2str(ECT2), ' para ', num2str(SNR_values(u)), ' dB']);

        P_senial = mean(c(1,:, 1).^2+c(2,:, 1).^2);
        SNR_est = 10 * log10(P_senial / P_noise(u));
        disp(['Relación señal-ruido (SNR) para ',num2str(SNR_values(u)),'dB de ruido ', num2str(SNR_est), ' dB']);
        
        desvio_std =median(abs(real(F(:))))/0.6745;
        Q_est = sqrt(2*gammaincinv(.997,1));
        espectrograma = abs(F).^2;
        umbral = Q_est * desvio_std;
        F_somb = F;
        F_somb (espectrograma < umbral) = 0;
        
        Plot_STFT(F_somb, t, f);
        hold on;
        plot(t, c(:,1, 1), 'r');
        plot(t, c(:,2, 1), 'b');
        legend('Cresta 1','Cresta 2');
        title(['Detección de cresta + umbralado con ' num2str(SNR_values(u)) 'dB de ruido']);
        hold off;
    
    end
    
    

%De esta forma siempre extraemos el primer modo





