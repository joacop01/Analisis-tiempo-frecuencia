function [STFT] = STFT_Gauss_diff(s, t, k)
    N = length(s);
    f = 0: 1: N-1;
    STFT = zeros(N/2, N); 
    for n = 1:N                          % Desplazamiento de la ventana
        g = -2*k*(t - t(n)).*exp(-k*(t - t(n)).^2);                 % Ventana gaussiana
        x_v = s .* g;                           % Multiplicación de la señal por la ventana
        x_v_f = fft(x_v).*exp(1i*2*pi*t(n).*f);
        x_v_f = x_v_f(1:N/2);
        STFT(:, n) = x_v_f; 
    end
end