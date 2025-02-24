function [x_estimado, v] = Wave_Shape(y, D, N, x, media)
    modulo = abs(y);
    fase = angle(y)/(2*pi);
    % Construir la matriz C
    C = zeros(2*D, N);           % Inicialización de C                   % Término constante (c_0)
    for d = 1:D
        C(d, :) = modulo.*cos(2 * pi *d * fase);    % Componentes coseno (c_d)
        C(d+D, :) = modulo.*sin(2 * pi * d * fase); % Componentes seno (d_d)
    end

    v = (x * C')*inv(C*C');       % v = (xC^T)(CC^T)^(-1)

    x_estimado = v * C;
    
    x_estimado = x_estimado + media;
end