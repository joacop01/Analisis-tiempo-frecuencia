function [env_max, env_min] = encontrar_envolventes(signal)
    N = length(signal);
    t = 1: length(signal);
    [pks_max, locs_max] = findpeaks(signal);
    
    % Encontrar mínimos locales
    [pks_min, locs_min] = findpeaks(-signal);
    pks_min = -pks_min; % Restaurar el signo para los mínimos

    % Ajustar bordes: agregar puntos inicial y final
    valor_borde = 0 ; % Usar promedio de la señal como valor de borde
    valor_max = mean(signal);
    locs_max = [1, locs_max, length(signal)];
    pks_max = [mean(signal), pks_max, mean(signal)];
    
    locs_min = [1, locs_min, length(signal)];
    pks_min = [mean(signal), pks_min, mean(signal)];
    
    % Interpolar las envolventes con spline
    env_max = interp1(locs_max, pks_max, t, 'spline');
    env_min = interp1(locs_min, pks_min, t, 'spline');
    
%     for k = 1: N
%         if(env_max(k) < signal(k))
%             env_max(k) = signal(k);
%         end
%         if(env_min(k) > signal(k))
%             env_min(k) = signal(k);
%         end
%     end
end