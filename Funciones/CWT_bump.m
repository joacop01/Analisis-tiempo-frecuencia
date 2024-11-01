function [CWT] = CWT_bump(s, f, sigma, scales, u)
    x_f = fftshift(fft(s));
    CWT = zeros(length(scales),length(s));
    for k = 1:length(scales)

        
        denominator = 1 - ((f*scales(k) - u).^2)/(sigma^2);  % Calcula el denominador de la función bump
        denominator(denominator == 0) = eps; % Asegura que no sea cero o negativo
        bump_analitica = exp(1 - (1 ./ (denominator))); % Función bump
        
        % Crear una envolvente en función del tiempo
        chi = (f*scales(k) >= u - sigma) & (f*scales(k) <= u + sigma); % Define la ventana temporal

        % Aplicar la envolvente
        bump_analitica = bump_analitica .* chi; % Aplicar la envolvente a la función bump
        bump_analitica(f*scales(k) < 0) = 0;  % Poner 0 donde f es negativo
        bump_analitica(f*scales(k) > 0) = 2*bump_analitica(f*scales(k) > 0);


        X_wav_f_analitica = x_f .* conj(bump_analitica);

        % Inversa de Fourier para obtener la CWT en el dominio del tiempo
        CWT(k, :) = ifft(X_wav_f_analitica);

            
    end
end