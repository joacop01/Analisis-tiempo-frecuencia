function [T] = Synchro_STFT(F, F_g, t, N)
    f = 0: 1 : N/2 -1;
    f = transpose(f);
    f_2 = repmat(f, 1, N);
    
    f_inst = f_2-((1/(2*pi)).*imag(F_g./ F));  % Frecuencia instantánea
    
    T = zeros(size(F));
    
    for n = 1: length(t)
        for k = 1:length(f)
            
            % Sincronización: reasignación de la energía
            a = f_inst(:, n);
            a = round(a);

            % Verificar si 'a' es un índice válido

            if a(k) >= 1 && a(k) <= length(f)
                if a(k) == f(k)
                    T(a(k), n) = T(a(k), n) +  (F(k, n));
                end
            end
        end

    end  
end
