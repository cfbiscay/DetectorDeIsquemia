clear;

derivacion_analizar = 'MLIII';

beta = 0.001;

load(['<file_path>' '\' derivacion_analizar]);  % <file_path> debe ser cambiado por la dirección donde se guardó qué registros y en qué derivación se utilizó la derivación MLIII (variables: deriv_registro, deriv_derivaciones)

% Recorro cada registro
for r = 36:length(deriv_registro(:,1))
    
    reg = deriv_registro(r,:); % nombre del registro
    d = deriv_derivaciones(r); % derivación de ese registro a analizar
    
    disp(reg);
    
    load(['<file_path>' '\' reg '.mat']); % <file_path> debe ser cambiado por la dirección donde se guardaron los archivos .mat generados a partir de las anotaciones de la base de datos (variable: anotaciones)  
    load(['<file_path>' '\' reg '.mat']); % <file_path> debe ser cambiado por la dirección donde se guardaron los archivos .mat generados a partir de los .ari (variable: marcas)
    load(['<file_path>' '\' reg '.mat']); % <file_path> debe ser cambiado por la dirección donde se guardaron las series con los parámetros extraídos
    
%% Genero un flag de anotaciones de la base de datos

    if ~isempty(anotaciones.inicio_fin)
        inicio_EP_real = anotaciones.x_derivacion.inicio_fin(:,1,d)*250;
        fin_EP_real = anotaciones.x_derivacion.inicio_fin(:,2,d)*250;
        extrema_EP_real = anotaciones.x_derivacion.extrema(:,1,d)*250;

        inicio_EP_real_2deriv = anotaciones.inicio_fin(:,1)*250;
        fin_EP_real_2deriv = anotaciones.inicio_fin(:,2)*250;
        extrema_EP_real_2deriv = anotaciones.extrema(:,1)*250;
    else
        inicio_EP_real1 = [];
        fin_EP_real1 = [];
        extrema_EP_real1 = [];

        inicio_EP_real_2deriv = [];
        fin_EP_real_2deriv = [];
        extrema_EP_real_2deriv = [];
    end
    
    latidos_eliminados = or(eliminar_latidos1,eliminar_latidos2);
    
    EP_flag_real = zeros(length(t),1);
    for i=1:length(inicio_EP_real)
        if fin_EP_real(i)==0
            break;
        end
        [~,aux1] = min(abs(t-inicio_EP_real(i)));
        [~,aux2] = min(abs(t-fin_EP_real(i)));
        EP_flag_real(aux1:aux2)=1;
    end

%% Filtro y normalizo las series de los parámetros y genero una única serie a partir de las anteriores

    Hil_filt22(1:6) = serie_Hilbert_sin_rechazo(1:6);
    for i = 7 : length(serie_ActHjorth_sin_rechazo)
        if latidos_eliminados(i+5)
            Hil_filt22(i-5:i) = Hil_filt22(i-6);
        else
            Hil_filt22(i) = serie_Hilbert_sin_rechazo(i);
        end
    end
    
    t_aux = removerows(t_sin_rechazo','ind',find(eliminar_latidos));
    
    [Hil_filt2,t22] = filtrado_y_normalizado(removerows(Hil_filt22','ind',find(eliminar_latidos)),t_aux,100);
    
    Act_filt22(1:6) = serie_ActHjorth_sin_rechazo(1:6);
    for i = 7 : length(serie_ActHjorth_sin_rechazo)
        if latidos_eliminados(i+5)
            Act_filt22(i-5:i) = Act_filt22(i-6);
        else
            Act_filt22(i) = serie_ActHjorth_sin_rechazo(i);
        end
    end
    
    t_aux = removerows(t_sin_rechazo','ind',find(eliminar_latidos));
    
    [Act_filt2,t2] = filtrado_y_normalizado(removerows(Act_filt22','ind',find(eliminar_latidos)),t_aux,100);
    
    Serie_param = mean([abs(Act_filt2) abs(Hil_filt2)]')'; % Serie que voy a usar para la detección
    
%% Busco el valor de nu para el umbral ajustando la serie a dos gaussianas

    GMModel = fitgmdist(Serie_param,2);
    idx = cluster(GMModel,Serie_param);
    cluster1 = Serie_param(idx == 1,:);
    cluster2 = Serie_param(idx == 2,:);
    gcluster1 = fitdist(cluster1,'Normal');
    gcluster2 = fitdist(cluster2,'Normal');
    x_values = 0:.01:10000;
    y1 = pdf(gcluster1,x_values);
    y2 = pdf(gcluster2,x_values);
    [max_y1,x_max_y1] = max(y1);
    [max_y2,x_max_y2] = max(y2);
    if x_max_y1 < x_max_y2
        x_start = x_max_y1;
        x_end = x_max_y2;
    else
        x_start = x_max_y2;
        x_end = x_max_y1;
    end
    [~,cruce] = min(abs(y1(x_start:x_end)-y2(x_start:x_end)));
    nu = x_values(cruce+x_start-1) - x_values(x_start);
    
%% Realizo la detección

    [EP_flag2,inicio_EP2,fin_EP2,baseline2,umbral2] = deteccion_EP(abs(Serie_param),t,beta,nu,inicio_EP_real,fin_EP_real,extrema_EP_real);

%% Grafico los resultados
    figure
    plot(t2,Serie_param)
    
    hold on
    plot(t2,baseline2)
    plot(t2,baseline2+nu)
    plot(t2,EP_flag2/2-1)
    plot(t2,EP_flag_real/3-1)
    hold off
  
%%
    clearvars Hil_filt22 Act_filt22
end