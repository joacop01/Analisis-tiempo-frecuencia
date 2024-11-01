addpath('../Funciones');
%Mexhat
% Parámetros de muestreo
Fs = 1000;             % Frecuencia de muestreo
Ts = 1/Fs;             % Periodo de muestreo
t = 0:Ts:1-Ts;             % Vector de tiempo
N = length(t);
f = -Fs/2:Fs/N: Fs/2-Fs/N; %Vector de frecuencias 

% Señal de prueba
% x = cos(2*pi*100*t);     % Ejemplo de señal senoidal
x1 = cos(2*pi*25*t + 2*pi*100*t.^2);
x2 = cos(2*pi*150*t + 2*pi*100*t.^2);
% 
x = x1 + x2;

% Parámetros de la wavelet sombrero mexicano
sigma = 1000;           % Parámetro de escala de la Gaussiana

% FFT de la señal
x_f = fftshift(fft(x));

% Definir las escalas (puedes ajustar los valores de escala)
scales = linspace(0,1,N);
% Inicializar la matriz de CWT
CWT = zeros(length(scales),length(x));

% Realizar la CWT en el dominio de la frecuencia

for k = 1:length(scales)
%     scale = scales(k);
    % Redefinir la variable de tiempo para la escala actual
%     t_scaled = (t-t(k)) / scale;
%     f_scaled = f*scale;
    % Definir la wavelet sombrero mexicano escalada
%     gauss = exp(-sigma*(t_scaled.^2)); 
%     mex = (-2*sigma*(2*sigma*t_scaled.^2-1) .* gauss)/(2*sigma);

    %mex_f = fftshift(fft(mex)); 
    mex_f = (2*pi^2*(f*scales(k)).^2)*sqrt(pi/sigma).*exp((-pi^2*(f*scales(k)).^2)/(sigma));
    % Multiplicación en frecuencia (señal * conj(wavelet escalada))
    X_wav_f = x_f .* conj(mex_f);

    % Inversa de Fourier para obtener la CWT en el dominio del tiempo
    CWT(k, :) = ifft(X_wav_f);
end

% Graficar la CWT


% Inicializar la matriz de CWT
CWT_analitica = zeros(length(scales),length(x));

% Realizar la CWT en el dominio de la frecuencia

for k = 1:length(scales)
    %length(scales)

    % Redefinir la variable de tiempo para la escala actual

    % Definir la wavelet sombrero mexicano escalada
    %gauss = exp(-sigma*t_scaled.^2); 
    %mex_analitica = (-2*sigma*(2*sigma*t_scaled.^2-1) .* gauss)/(2*sigma);
    % FFT de la wavelet escalada
    %Mex_f_analitica = fftshift(fft(mex_analitica));
    Mex_f_analitica = (2*pi^2*(f*scales(k)).^2)*sqrt(pi/sigma).*exp((-pi^2*(f*scales(k)).^2)/(sigma));
    Mex_f_analitica(f < 0) = 0;  % Poner 0 donde f es negativo
    Mex_f_analitica(f > 0) = 2*Mex_f_analitica(f > 0);
    X_wav_f_analitica = x_f .* conj(Mex_f_analitica);

    % Inversa de Fourier para obtener la CWT en el dominio del tiempo
    CWT_analitica(k, :) = ifft(X_wav_f_analitica);
end

%Graficar la CWT
figure;
subplot(211);
imagesc(t, linspace(0,N),abs(CWT_analitica));  % Graficar la magnitud de la CWT
xlabel('Tiempo');
ylabel('Escalas');
title('CWT con ondita mexhat (analitica)');
colorbar;
subplot(212);
imagesc(t, linspace(0, N), abs(CWT));  % Graficar la magnitud de la CWT
xlabel('Tiempo');
ylabel('Escalas');
title('CWT con ondita mexhat (real)');
colorbar;
