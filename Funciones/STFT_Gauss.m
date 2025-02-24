function [STFT] = STFT_Gauss(s, t, k)
    N = length(s);

    f = 0: 1 : N-1;
    STFT = zeros(round(N/2), N); 
    for n = 1:N                          
        g = exp(-k*(t - t(n)).^2);                 
        x_v = s .* g;                           
        x_v_f = fft(x_v).*exp(1i*2*pi*t(n).*f);
        x_v_f = x_v_f(1:N/2);
        STFT(:, n) = x_v_f; 
    end
end