function [CWT] = CWT_morlet_diff(s, f,  sigma, scales, f0)
    x_f = fftshift(fft(s));
    N = length(s);
    CWT = zeros(length(scales),length(s));
    for k = 1:length(scales)

        
        morlet_f_analitica = sqrt(pi/sigma)*exp((-pi^2*(f*scales(k)-f0).^2)/sigma);
        morlet_f_analitica = ((-2*pi^2*(f*scales(k)-f0))/(sigma)).*morlet_f_analitica;
        morlet_f_analitica(f*scales(k) < 0) = 0;  % Poner 0 donde f es negativo
        morlet_f_analitica(f*scales(k) > 0) = 2*morlet_f_analitica(f*scales(k) > 0);
        X_wav_f_analitica = x_f .* conj(morlet_f_analitica);

        % Inversa de Fourier para obtener la CWT en el dominio del tiempo
        CWT(k, :) = ifft(X_wav_f_analitica);
    end
end