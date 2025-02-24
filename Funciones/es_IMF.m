function is_imf = es_IMF(signal)
    % Verifica si una señal es una IMF según las condiciones:
    % 1. Las cantidades de extremos y cruces por cero difieren a lo sumo en 1.
    % 2. El valor medio de las envolventes superior e inferior es cercano a 0.
    
    % Condición 1: Cantidad de extremos y cruces por cero
    num_extremos = sum(diff(sign(diff(signal))) ~= 0);
    num_cruces_cero = sum(diff(sign(signal)) ~= 0);
    
    % Condición 2: Valor medio de las envolventes
    [env_max, env_min] = encontrar_envolventes(signal);
    mean_env = (env_max + env_min) / 2;
    mean_close_to_zero = mean((abs(mean_env))) < 0.01;
    
    % Validar ambas condiciones
    is_imf = abs(num_extremos - num_cruces_cero) <= 1 && mean_close_to_zero;
end