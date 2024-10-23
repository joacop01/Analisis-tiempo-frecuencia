function [STFT] = STFT_Hann(s, t, T)
    N = length(s);
    STFT = zeros(N, N); 
    for n = 1:N
        
        chi = (t - t(n) >= -T/2) & (t - t(n) <= T/2);
        g = (1 + cos(2*pi*(t-t(n))/T)) .* chi;
        x_v = s .* g;                           % Multiplicación de la señal por la ventana
        x_v_f = fftshift(fft(x_v)).*exp(1i*2*pi*t/N);
        %x_v_f = x_v_f(N/2+1:end);
        STFT(:, n) = x_v_f;                     % FFT y fftshift para centrar la frecuencia
    end
end