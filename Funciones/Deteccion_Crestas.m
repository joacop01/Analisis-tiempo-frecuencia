function [crestas, F_aux] = Deteccion_Crestas(F, n, N, cant_crestas, Q)
   
    indices = round(linspace(0.5/(n+1), 1 - 0.5/(n+1), n) * N);
    c = zeros(N, cant_crestas, length(indices));
    energia = zeros(cant_crestas, 1);
    indice_max = zeros(cant_crestas, 1); 
    crestas = zeros(N, cant_crestas);
    for i = 1: length(indices)
        F_aux = F;
        for u = 1: cant_crestas
            [~, c(indices(i), u, i)] = max(abs(F_aux(:,indices(i)).^2));

            for a = indices(i)+1:N
                I_a = c(a-1, u, i) - Q: c(a-1, u, i) + Q;

                I_a = max(I_a(1), 1):min(I_a(end), N/2);

                [~, indice_local_a] = max(abs(F_aux(I_a, a)).^2);

                c(a, u, i) =  I_a(indice_local_a);

            end

            for b = indices(i)-1:-1:1
                I_b = c(b+1, u, i) - Q: c(b+1, u, i) + Q;

                I_b = max(I_b(1), 1):min(I_b(end), N/2);

                [~, indice_local_b] = max(abs(F_aux(I_b, b)).^2);

                c(b, u, i) = I_b(indice_local_b);
            end

            %Borramos la cresta recien detectada
            for k = 1:N
                I = c(k, u, i) - Q : c(k, u, i) + Q;
                I = max(I(1), 1):min(I(end), N/2);
                F_aux(I, k) = 0;
            end
            energia_aux = sum(abs(F(c(:, u, i)))).^2;
            energia_borrada = sum(abs(F_aux(c(:, u, i))));
            if energia_borrada < 1000
                if energia_aux > energia(u)
                    energia(u) = energia_aux;
                    indice_max(u) = i;
                end
            else
                indice_max(u) = 1;
            end
        end

    end

    for w = 1:cant_crestas
        crestas(:, w) = c(:, w, indice_max(w));
    end
    
end
