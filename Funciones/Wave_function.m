function [wave_func, t] = Wave_function(v, D, N, Fs)
    Ts = 1/ Fs;
    t = 0: Ts: 1 - Ts;
    C = zeros (2*D, N);
    for d = 1: D
        C(d, :) = cos(2 * pi *d * t);    % Componentes coseno (c_d)
        C(D+d, :) = sin(2 * pi * d * t); % Componentes seno (d_d)
    end
    
    wave_func = v * C;
end