function [STFT] = STFT_Gauss(s, t, k)
    N = length(s);
%     f = 0 : 1/N : 1 - 1/N;
    f = 0: 1 : N-1;
    STFT = zeros(N/2, N); 
    for n = 1:N                          % Desplazamiento de la ventana
        g = exp(-k*(t - t(n)).^2);                 % Ventana gaussiana
        x_v = s .* g;                           % Multiplicación de la señal por la ventana
        %x_v_f = fftshift(fft(x_v)).*exp(1i*2*pi*t/N);
        x_v_f = fft(x_v).*exp(1i*2*pi*t(n).*f);
        x_v_f = x_v_f(1:N/2);
        STFT(:, n) = x_v_f; 
    end
end