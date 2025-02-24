function [IMFs] = EMD(signal, modos)
    addpath('../Funciones');

    r = signal;
    IMFs = zeros(modos, length(signal));
    j = 1;
    k = 1;
    d = []; % Matriz para almacenar todas las funciones descomponedoras
    suma_d = zeros(1, length(signal)); % Inicializa el acumulador de todas las d

    while(k < modos)
        while (true)  
            [env_max, env_min] = encontrar_envolventes(r);

            m = (env_max + env_min) / 2;

            d(j, :) = r - m;

            if(es_IMF(d(j,:)))
                suma_d = suma_d + d(j,:); % Acumular la IMF
                IMFs(k, :) = d(j, :);
                a = signal - suma_d;
                r = a;
                k = k +1;
                j = 1;
                break;
            else
                r = d(j, :);
                j = j + 1;
            end
        end
        IMFs(modos, :) = r;

end