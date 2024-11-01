function [CWT] = CWT_mexhat(s, f, sigma, scales)
    x_f = fftshift(fft(s));
    CWT = zeros(length(scales),length(s));
    for k = 1:length(scales)
        

        Mex_f_analitica = (2*pi^2*(f*scales(k)).^2)*sqrt(pi/sigma).*exp((-pi^2*(f*scales(k)).^2)/(sigma));
        Mex_f_analitica(f*scales(k) < 0) = 0;  % Poner 0 donde f es negativo
        Mex_f_analitica(f*scales(k) > 0) = 2*Mex_f_analitica(f*scales(k) > 0);

        X_wav_f_analitica = x_f .* conj(Mex_f_analitica);


        % Inversa de Fourier para obtener la CWT en el dominio del tiempo
        CWT(k, :) = ifft(X_wav_f_analitica);
        
    end
end