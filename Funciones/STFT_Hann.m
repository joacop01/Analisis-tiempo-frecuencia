function [STFT] = STFT_Hann(s, t, T)
    N = length(s);
    STFT = zeros(round(N/2), N); 
    f = 0: 1 : N-1;
    for n = 1:N
        chi = (t - t(n) >= -T/2) & (t - t(n) <= T/2);
        g = (1 + cos(2*pi*(t-t(n))/T)) .* chi;
        x_v = s .* g;                          
        x_v_f = fft(x_v).*exp(1i*2*pi*t(n).*f);
        x_v_f = x_v_f(1:N/2);
        STFT(:, n) = x_v_f;                    
    end
end