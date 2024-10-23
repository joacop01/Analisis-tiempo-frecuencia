Fs = 1000;                  % Frecuencia de muestreo
Ts = 1/Fs;                  % Intervalo de muestreo
t = 0:Ts:1-Ts;              % Vector de tiempo
N = length(t);              % Longitud de la señal
k = 200;                    % Parámetro de la ventana gaussiana
len = 0.1:0.01:1;           % Diferentes longitudes de la ventana
sigma = 18 ./ len.^2;       % Parámetro sigma para cada longitud de ventana
f1 = 100;                   % Frecuencia de la primera señal
f2 = 150;                   % Frecuencia de la segunda señal
alpha1 = 2;                 % Parámetro alpha para la entropía de Rényi (caso 1)
alpha2 = 3;                 % Parámetro alpha para la entropía de Rényi (caso 2)

% x1 = cos(2*pi*f1*t);                  % Primera señal
% x2 = cos(2*pi*f2*t);                  % Segunda señal
x1 = cos(2*pi*100*t + 2*pi*100*t.^2);
x2 = cos(2*pi*150*t + 2*pi*100*t.^2);
x = x1 + x2;                          % Señal combinada

P_signal = mean(x.^2);                % Potencia de la señal original
% Escalar el ruido para diferentes SNR
% Escalar el ruido para diferentes SNR
SNR_values = [0, -5, -10];            % Valores de SNR en dB
ruido_gaussiano = zeros(size(x));     % Inicializar el ruido
x_ruido = zeros(length(SNR_values), length(x));  % Inicializar matriz para las señales ruidosas

for idx = 1:length(SNR_values)
    snr_db = SNR_values(idx);          % Obtener el valor de SNR actual
    
    % Calcular la potencia del ruido deseada
    P_noise = P_signal / 10^(snr_db / 10);
    
    % Generar ruido gaussiano con la potencia deseada
    ruido_gaussiano = sqrt(P_noise) * randn(size(x));  % Escalar el ruido
    
    % Añadir ruido a la señal y almacenar en la columna correspondiente
    x_ruido(idx, :) = x + ruido_gaussiano;  % Transponer para que cada columna corresponda a una señal ruidosa
    P_noise = mean(ruido_gaussiano.^2);  % Potencia del ruido
    SNR = 10 * log10(P_signal / P_noise);

    % Mostrar el valor de la SNR (debería ser cercano a 0 dB)
    disp(['Relación señal-ruido (SNR): ', num2str(SNR), ' dB']);
end


% Parámetros para la STFT
STFT = zeros(N, N);                    % Matriz para almacenar el resultado de la STFT

% Inicialización de los vectores para almacenar la entropía de Rényi
renyi_alpha1 = zeros(1, length(len)); 
renyi_alpha2 = zeros(1, length(len)); 
for u = 1:length(SNR_values)
    for l = 1:length(len)
        for n = 1:N
            a = n / Fs;                           % Desplazamiento de la ventana
            g = exp(-sigma(l) * (t - a).^2);      % Ventana gaussiana
            x_v = x_ruido(u,:);
            x_v = x_v .* g;                         % Multiplicación de la señal por la ventana
            x_v_f = fftshift(fft(x_v)).*exp(1i*2*pi*t/N);           % FFT de la señal segmentada
            STFT(:, n) = x_v_f;                   % Almacenar en la matriz de STFT
        end

        % Cálculo del espectrograma
        spec = abs(STFT).^2;                      % Potencia espectral (magnitud al cuadrado)

        % Normalización del espectrograma
        spec_norm = spec ./ sum(sum(spec));         % Normalización para que el espectrograma sea una distribución de probabilidad

        % Cálculo de la entropía de Rényi para alpha = 2
        renyi_alpha1(l) = (1 / (1 - alpha1)) * log2(sum(sum(spec_norm.^alpha1)));

        % Cálculo de la entropía de Rényi para alpha = 3
        renyi_alpha2(l) = (1 / (1 - alpha2)) * log2(sum(sum(spec_norm.^alpha2)));

    end
    
    figure;
    % Subplot para la señal original
    subplot(2, 1, 1);
    plot(t, x, 'LineWidth', 1.5);
    title('Señal Original');
    xlabel('Tiempo (s)');
    ylabel('Amplitud');
    grid on;

    % Subplot para la señal ruidosa con SNR = 0 dB
    subplot(2, 1, 2);
    plot(t, x_ruido(u, :), 'LineWidth', 1.5);
    title(['Señal Ruidosa (SNR =',num2str(SNR_values(u)),'dB)']);
%     title(['Espectrograma de la señal ruidosa (SNR = ', num2str(SNR_values(u)), ' dB) - Longitud Óptima de Ventana']);
    xlabel('Tiempo (s)');
    ylabel('Amplitud');
    grid on;

    % % Subplot para la señal ruidosa con SNR = -5 dB
    % subplot(4, 1, 3);
    % plot(t, x_ruido(2, :), 'LineWidth', 1.5);
    % title('Señal Ruidosa (SNR = -5 dB)');
    % xlabel('Tiempo (s)');
    % ylabel('Amplitud');
    % grid on;
    % 
    % % Subplot para la señal ruidosa con SNR = -10 dB
    % subplot(4, 1, 4);
    % plot(t, x_ruido(3, :), 'LineWidth', 1.5);
    % title('Señal Ruidosa (SNR = -10 dB)');
    % xlabel('Tiempo (s)');
    % ylabel('Amplitud');
    % grid on;

    % Graficar la entropía de Rényi en función de la longitud de la ventana
    figure;
    plot(len, renyi_alpha1, '-o', 'DisplayName', ['α = ', num2str(alpha1)]);
    hold on;
    plot(len, renyi_alpha2, '-s', 'DisplayName', ['α = ', num2str(alpha2)]);
    xlabel('Longitud de la ventana');
    ylabel('Entropía de Rényi');
    title('Entropía de Rényi vs Longitud de la ventana');
    legend show;

    % Encontrar el mínimo global para alpha = 2 y alpha = 3
    [min_renyi1, index_min_renyi1] = min(renyi_alpha1);
    [min_renyi2, index_min_renyi2] = min(renyi_alpha2);

    % Mostrar resultados del mínimo
    disp(['Mínimo global de la entropía de Rényi (α = ', num2str(alpha1), '): ', num2str(min_renyi1), ' en longitud de ventana: ', num2str(len(index_min_renyi1)), 's']);
    disp(['Mínimo global de la entropía de Rényi (α = ', num2str(alpha2), '): ', num2str(min_renyi2), ' en longitud de ventana: ', num2str(len(index_min_renyi2)),'s']);
    longitud_ven = index_min_renyi2;
    sigma_op = 18/longitud_ven^2;
    for n = 1:N
        a = n / Fs;                           % Desplazamiento de la ventana
        g = exp(-sigma_op* (t - a).^2);      % Ventana gaussiana
        x_v = x_ruido(u,:) .* g;                         % Multiplicación de la señal por la ventana
        x_v_f = fftshift(fft(x_v));           % FFT de la señal segmentada
        STFT(:, n) = x_v_f;                   % Almacenar en la matriz de STFT
    end
    k = 0:Fs/N: Fs/2 - Fs/N;                    % Frecuencias positivas

    % Graficar la señal original y su STFT
    figure;
    subplot(2,1,1);
    plot(t, real(x_ruido(u,:)));
    title('Señal en el tiempo');
    xlabel('Tiempo (s)');
    ylabel('Amplitud');
    %hold on;
    %plot(t, g);
    %title('Ventana gaussiana');
    %hold off;

    % Mostrar el espectrograma (STFT en función del tiempo y frecuencia)
    subplot(2,1,2);
    imagesc(t, k, abs(STFT(N/2+1:end,:)));   % Magnitud de la STFT

    % Invertir el eje de frecuencias (frecuencias más altas abajo)
    set(gca, 'YDir', 'reverse');  % Ahora las frecuencias altas estarán abajo

    xlabel('Tiempo (s)');
    ylabel('Frecuencia (Hz)');
    title(['Espectrograma de la señal ruidosa (SNR = ', num2str(SNR_values(u)), ' dB) - Longitud Óptima de Ventana']);
    colorbar;
    
        % Marginalización de la STFT respecto al tiempo
    % Calcular la potencia espectral (espectrograma)
    espectrograma = abs(STFT(N/2+1:end,:)).^2; % Potencia espectral

    % Sumar la potencia espectral a lo largo del tiempo
    marginal = sum(espectrograma, 2); % Sumar a lo largo del tiempo (columnas)

    % Normalizar la marginal
    marginal_normalizada = marginal / sum(marginal); % Normalización para que sea una distribución
    
      % Calcular la FFT de la señal ruidosa
    fft_s = fftshift(fft(x_ruido(u,:)));
    fft_s_norm = abs(fft_s).^2 / sum(abs(fft_s).^2); % Normalizar la FFT

    % Graficar la marginalización
    figure;
    plot(k, marginal_normalizada, 'LineWidth', 1.5, 'DisplayName', 'Marginalización STFT');
    hold on;
    plot(k, fft_s_norm(N/2+1:end), 'LineWidth', 1.5, 'DisplayName', 'FFT de la Señal Ruidosa');
    
    xlabel('Frecuencia (Hz)');
    ylabel('Densidad Espectral de Potencia Normalizada');
    title(['Marginalización de la STFT respecto al tiempo (SNR = ', num2str(SNR_values(u)), ' dB)']);
    grid on;

    % Agregar leyenda para identificar las curvas
    legend show;

    hold off;
    
end