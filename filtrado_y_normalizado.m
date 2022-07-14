function [serie_filt_y_norm,t_nuevo] = filtrado_y_normalizado(serie,t,n)
% Filtra, normaliza y remuestrea la serie

%% entrada:
%   serie: serie que quiero normalizar y filtrar
%   t: vector de tiempo
%   n: muestras


%% salida
%   serie_filt_y_norm: serie remuestreada, normalizada y filtrar
%   t_nuevo: vector de tiempo de la serie remuestreada

%% 
%%

%% Remuestreo la serie
[serie_filtrada,t_nuevo] = resample(serie,t,1/250);   

%% Filtrado de la serie
serie_filtrada = smooth(serie_filtrada,n);
   
%% Normalizado de la serie (uso los primeros 100 latidos para normalizar)
serie_basal = serie_filtrada(1:100);
GMModel = fitgmdist(serie_basal,1);
idx = cluster(GMModel,serie_basal);
cluster1 = serie_basal(idx == 1);
gcluster1 = fitdist(cluster1,'Normal');
if max(serie_basal)<0.01
    x_values = min(serie_basal)*2:.01:max(serie_basal);
else
    x_values = 0:.01:max(serie_basal)*2;
end

y1 = pdf(gcluster1,x_values);

[a,b] = max(y1);        
c = b - 1 + find(y1(b:end)<a*0.005,1);
if isempty(c)
    c=length(y1);
end

serie_filt_y_norm = (serie_filtrada - x_values(b)) / x_values(c);