function [cos_sint, v] = Sintetic_cos_func(a, b, D, modulo, t)
    N = length(t);
    C = zeros((2*D), N);
    fase = a*t + (b/(2*pi))*sin(2*pi*t);
    for d = 1:D
        C(d, :) = modulo*cos(2 * pi *d * fase);    % Componentes coseno (c_d)
        C(D + d, :) = modulo*sin(2 * pi * d * fase); % Componentes seno (d_d)
    end
    v = rand(1, 2*D);
    v(1) = 1;

    cos_sint = v * C;
end