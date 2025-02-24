function [CWT] = CWT_mexhat(s, f, sigma, scales)
    x_f = fftshift(fft(s));
    CWT = zeros(length(scales),length(s));
    figure;
for k = 1:length(scales)

    f_scaled = f* scales(k);
    Mex_f_analitica = 2*pi^2*(f_scaled.^2)*sqrt(pi/sigma).*exp(-pi^2*(f_scaled.^2)/(sigma));
    Mex_f_analitica(f < 0) = 0;  % Poner 0 donde f es negativo
    Mex_f_analitica(f > 0) = 2*Mex_f_analitica(f > 0);
    X_wav_f_analitica = x_f .* conj(Mex_f_analitica);
    plot(f, Mex_f_analitica);
    hold on;

    % Inversa de Fourier para obtener la CWT en el dominio del tiempo
    CWT(k, :) = ifft(X_wav_f_analitica);
end
    hold off;
end