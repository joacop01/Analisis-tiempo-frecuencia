function [y] = Reconstruccion_Cresta(c, b, N, F, cant_crestas)
   y = zeros(1, N); 
   for k = 1: cant_crestas
       for i =1:N

        I = c(i, k) - b : c(i, k) + b;

        I = max(I(1), 1):min(I(end), N/2);

        y(i) = y(i) + sum(F(I,i));
       end
    end
end