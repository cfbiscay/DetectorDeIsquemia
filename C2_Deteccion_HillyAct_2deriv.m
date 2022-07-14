clear

record = ['e0103','e0104','e0105','e0106','e0107','e0108','e0110', ...
    'e0111','e0112','e0113','e0114','e0115','e0116','e0118','e0119', ...
    'e0121','e0122','e0123','e0124','e0125','e0126','e0127','e0129', ...
    'e0133','e0136','e0139','e0147','e0148','e0151','e0154','e0155', ...
    'e0159','e0161','e0162','e0163','e0166','e0170','e0202','e0203', ...
    'e0204','e0205','e0206','e0207','e0208','e0210','e0211','e0212', ...
    'e0213','e0302','e0303','e0304','e0305','e0306','e0403','e0404', ...
    'e0405','e0406','e0408','e0409','e0410','e0411','e0413','e0415', ...
    'e0417','e0418','e0501','e0509','e0515','e0601','e0602','e0603', ...
    'e0604','e0605','e0606','e0607','e0609','e0610','e0611','e0612', ...
    'e0613','e0614','e0615','e0704','e0801','e0808','e0817','e0818', ...
    'e1301','e1302','e1304']; 

% Registros a los que les voy a extraer los parámetros (desde primer_registro hasta ultimo_registro)
primer_registro = 'e0103'; 
ultimo_registro = 'e1304';
in = strfind(record,primer_registro);
fi = strfind(record,ultimo_registro);

beta = 0.001;

% Recorro cada registro
for r = in : 5 : fi
    reg = record(r:r+4); % nombre del registro
    
    load(['<file_path>' '\' reg '.mat']); % <file_path> debe ser cambiado por la dirección donde se guardaron los archivos .mat generados a partir de las anotaciones de la base de datos (variable: anotaciones)  
    load(['<file_path>' '\' reg '.mat']); % <file_path> debe ser cambiado por la dirección donde se guardaron los archivos .mat generados a partir de los .ari (variable: marcas)
    load(['<file_path>' '\' reg '.mat']); % <file_path> debe ser cambiado por la dirección donde se guardaron las series con los parámetros extraídos
    
    figure(1)
    
    % Recorro cada derivación de cada registro
    for d = 1 : 2
        
        disp(['registro:' reg ' - d: ' num2str(d) ' - ruido: ' num2str(ruido(d)) '% - r: ' num2str(r)]);
        
        if d==1
            t = t1;
            t_sin_rechazo = t_sin_rechazo1;
            serie_Hilbert_sin_rechazo = serie_Hilbert_sin_rechazo1;
            serie_ActHjorth_sin_rechazo = serie_ActHjorth_sin_rechazo1;
            eliminar_latidos = eliminar_latidos1;
            dd = 2;
        elseif d==2
            t = t2;
            t_sin_rechazo = t_sin_rechazo2;
            serie_Hilbert_sin_rechazo = serie_Hilbert_sin_rechazo2;
            serie_ActHjorth_sin_rechazo = serie_ActHjorth_sin_rechazo2;
            eliminar_latidos = eliminar_latidos2;
            dd = 1;
        end

%% Genero un flag de anotaciones de la base de datos
        if ~isempty(anotaciones.inicio_fin)
            inicio_EP_real = anotaciones.x_derivacion.inicio_fin(:,1,d)*250;
            fin_EP_real = anotaciones.x_derivacion.inicio_fin(:,2,d)*250;
            extrema_EP_real = anotaciones.x_derivacion.extrema(:,1,d)*250;
            
            inicio_EP_real_otraDeriv = anotaciones.x_derivacion.inicio_fin(:,1,dd)*250;
            fin_EP_real_otraDeriv = anotaciones.x_derivacion.inicio_fin(:,2,dd)*250;
            extrema_EP_real_otraDeriv = anotaciones.x_derivacion.extrema(:,1,dd)*250;

            inicio_EP_real_2deriv = anotaciones.inicio_fin(:,1)*250;
            fin_EP_real_2deriv = anotaciones.inicio_fin(:,2)*250;
            extrema_EP_real_2deriv = anotaciones.extrema(:,1)*250;
        else
            inicio_EP_real = [];
            fin_EP_real = [];
            extrema_EP_real = [];
            
            inicio_EP_real_otraDeriv = [];
            fin_EP_real_otraDeriv = [];
            extrema_EP_real_otraDeriv = [];

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
        
        EP_flag_real_otraDeriv = zeros(length(t),1);
        for i=1:length(inicio_EP_real_otraDeriv)
            if fin_EP_real_otraDeriv(i)==0
                break;
            end
            [~,aux1] = min(abs(t-inicio_EP_real_otraDeriv(i)));
            [~,aux2] = min(abs(t-fin_EP_real_otraDeriv(i)));
            EP_flag_real_otraDeriv(aux1:aux2)=1;
        end

%% Filtro y normalizo las series de los parámetros y genero una única serie a partir de las anteriores

        Hil_filt22(1:6) = serie_Hilbert_sin_rechazo(1:6);
        for i = 7 : length(serie_Hilbert_sin_rechazo)-5
            if latidos_eliminados(i+5)
                Hil_filt22(i-5:i) = Hil_filt22(i-6);
            else
                Hil_filt22(i) = serie_Hilbert_sin_rechazo(i);
            end
        end        
        Hil_filt22(i+1:i+5) = serie_Hilbert_sin_rechazo(i+1:i+5);

        t_aux = removerows(t_sin_rechazo,'ind',find(eliminar_latidos));
        
        [Hil_filt2,~] = filtrado_y_normalizado(removerows(Hil_filt22','ind',find(eliminar_latidos)),t_aux,100);

        Act_filt22(1:6) = serie_ActHjorth_sin_rechazo(1:6);
        for i = 7 : length(serie_ActHjorth_sin_rechazo)-5
            if latidos_eliminados(i+5)
                Act_filt22(i-5:i) = Act_filt22(i-6);
            else
                Act_filt22(i) = serie_ActHjorth_sin_rechazo(i);
            end
        end
        Act_filt22(i+1:i+5) = serie_ActHjorth_sin_rechazo(i+1:i+5);
        
        t_aux = removerows(t_sin_rechazo,'ind',find(eliminar_latidos));
        
        [Act_filt2,~] = filtrado_y_normalizado(removerows(Act_filt22','ind',find(eliminar_latidos)),t_aux,100);

        Param = mean([abs(Act_filt2) abs(Hil_filt2)]')';

        % Busco el valor de nu para el umbral
        GMModel = fitgmdist(Param,2);
        idx = cluster(GMModel,Param);
        cluster1 = Param(idx == 1,:);
        cluster2 = Param(idx == 2,:);
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

        [EP_flag2,inicio_EP2,fin_EP2,baseline2,umbral2] = deteccion_EP(abs(Param),t,beta,nu,inicio_EP_real,fin_EP_real,extrema_EP_real);

%% Grafico los resultados
        subplot(2,1,d)
        plot(t,Param)

        hold on
        plot(t,baseline2)
        plot(t,baseline2+nu)
        plot(t,EP_flag_real_otraDeriv/3-1,'black')
        plot(t,EP_flag_real/3-1,'g')
        plot(t,EP_flag2/2-1,'r')
        hold off

%%
        clearvars Hil_filt22 Act_filt22
    end
end