function [env_max, env_min] = calcular_envolventes_mod(signal, t)

    N = length(signal);
    % Detectar máximos y mínimos locales con corrección de bordes
    [pks_max, locs_max] = findpeaks(signal, t); % Máximos locales
    [pks_min, locs_min] = findpeaks(-signal, t); % Mínimos locales (invertir la señal)
    pks_min = -pks_min; % Revertir el signo de los mínimos
    
    % Agregar los extremos iniciales y finales si no están incluidos
    if t(1) < locs_max(1) || signal(1) > signal(2)
        locs_max = [t(1), locs_max];
        pks_max = [signal(1), pks_max];
    end
    if t(end) > locs_max(end) || signal(end) > signal(end-1)
        locs_max = [locs_max, t(end)];
        pks_max = [pks_max, signal(end)];
    end
    
    if t(1) < locs_min(1) || signal(1) < signal(2)
        locs_min = [t(1), locs_min];
        pks_min = [signal(1), pks_min];
    end
    if t(end) > locs_min(end) || signal(end) < signal(end-1)
        locs_min = [locs_min, t(end)];
        pks_min = [pks_min, signal(end)];
    end

    % Interpolar las envolventes
    env_max = interp1(locs_max, pks_max, t, 'pchip'); % Envolvente máxima
    env_min = interp1(locs_min, pks_min, t, 'pchip'); % Envolvente mínima
    
    for k = 1: N
        if(env_max(k) < signal(k))
            env_max(k) = signal(k);
        end
        if(env_min(k) > signal(k))
            env_min(k) = signal(k);
        end
    end
end