function [marcasQRS_post_rechazos, latidos_rechazados, latidos_a_rechazar_en_marcasQRS] = rechazo_por_baseline(t,baseline,marcasQRS_normales)
% Rechaza latidos en donde se produzca un salto abrupto de la línea de base

%% entrada:
%   t: vector de tiempo
%   baseline: linea de base del ECG
%   marcasQRS_normales: indica el QRS de los latidos considerados normales (según las marcas de ARISTOTLE)



%% salida
%   marcasQRS_post_rechazos: indica el QRS de los latidos considerados normales y con los latidos eliminados por cambios abruptos en la línea de base
%   latidos_rechazados: indican los latidos rechazados por cambios abruptos de la línea de base
%   latidos_a_rechazar_en_marcasQRS: la ubicación de los latidos que se rechazan por cambios abruptos de la línea de base

%% 
%%

baseline_med = zeros(length(marcasQRS_normales));
rechazo = zeros(size(marcasQRS_normales));

baseline_med(1) = mean(baseline(marcasQRS_normales(1,1):marcasQRS_normales(2,1)-1));

for k=2:length(marcasQRS_normales)-1
    baseline_med(k) = mean(baseline(marcasQRS_normales(k):marcasQRS_normales(k+1)-1));
    if(abs(baseline_med(k)-baseline_med(k-1))>20)
        rechazo(k-1:k+1) = 1;
    end
end                                                              

latidos_a_rechazar_en_marcasQRS = find(rechazo);
latidos_a_rechazar_en_t = zeros(size(latidos_a_rechazar_en_marcasQRS));

for i = 1 : length(latidos_a_rechazar_en_marcasQRS)
    [~,latidos_a_rechazar_en_t(i)] = min(abs(t-marcasQRS_normales(latidos_a_rechazar_en_marcasQRS(i))));
end

marcasQRS_post_rechazos = removerows(marcasQRS_normales,'ind',latidos_a_rechazar_en_marcasQRS);
latidos_rechazados = marcasQRS_normales(latidos_a_rechazar_en_marcasQRS);