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

direcGuardado = '<file_path>'; % <file_path> debe ser cambiado por la dirección donde voy a guardar los parámetros

modos_de_interes = [8 9 10]; % IMFs seleccionados de interés
ruido = [0 0];

% Recorro cada registro
for r = in : 5 : fi
    reg = record(r:r+4);
    disp(reg);
    
    load('<file_path>'); % <file_path> debe ser cambiado por la dirección donde se guardaron los registros de la base de datos (variable: reg)
    load('<file_path>'); % <file_path> debe ser cambiado por la dirección donde se guardaron los archivos generados a partir de los .ari (variable: marcas)
         
%% Rechazo los latidos anormales (según las marcas ARI)

    eliminar_latidos_QRS = zeros(size(marcas.time));

    for i=1:10
        if marcas.anntyp(i) ~= 'N'
            eliminar_latidos_QRS(1:i+10) = 1;
        end    
    end
    for i=11:length(marcas.time)-10
        if marcas.anntyp(i) ~= 'N'
            eliminar_latidos_QRS(i-10:i+10) = 1;
        end
    end
    for i=length(marcas.time)-10:length(marcas.time)
        if marcas.anntyp(i) ~= 'N'
            eliminar_latidos_QRS(i-10:end) = 1;
        end
    end

    marcasQRS_normales = removerows(marcas.time,find(eliminar_latidos_QRS));
    
    for d = 1 : 2
        
%% Rechazo por el ruido en la linea de base
    
        load('<file_path>'); % <file_path> debe ser cambiado por la dirección donde se guardaron marcas del delineador (variable: pos)
    
        baseline = obtener_baseline(registro.ECG(:,d+1),[pos.Poff]',[pos.QRSon]',1,length(registro.ECG(:,d)),0);
        t_muestras = 1:length(baseline(:,1));
        
        [marcasQRS_post_rechazos, latidos_rechazados, ~] = rechazo_por_baseline(t_muestras,baseline,marcasQRS_normales);
        
        eliminar_latidos = eliminar_latidos_QRS;
        
        for i = 1:length(latidos_rechazados)
            eliminar_latidos(find(latidos_rechazados(i)==marcas.time)) = 1;            
        end
        
        ruido(d) = sum(eliminar_latidos)/length(eliminar_latidos)*100;
   
%%  Obtengo las series de los parámetros: Actividad de Hjorth y Transformada de Hilbert

        load('<file_path>'); % <file_path> debe ser cambiado por la dirección donde se guardaron los resultados de EMD (variable: imf)
    
        rec_imf = sum(imf(modos_de_interes,:)); % Hago la sumatorea de los modos de interés
        
        [~, a] = FAhilbert(imf(modos_de_interes,:),1/registro.FS); %% Calculo la Transformada de los IMFs de interés
        serie_Hilbert_aux = sum(a); % Hago la sumatoria de la amplitud de la Transformada de Hilbert

        k=1;
        for i=6:1:length(marcas.time)-5
            serie_ActHjorth_sin_rechazo(k) = var(rec_imf(marcas.time(i-5):marcas.time(i+5)));          
            serie_Hilbert_sin_rechazo(k) = mean(serie_Hilbert_aux(marcas.time(i):marcas.time(i+1)));
            t_sin_rechazo(k) = marcas.time(i);
            k = k + 1;
        end
        
        eliminar_latidos_aux = eliminar_latidos(1:length(serie_ActHjorth_sin_rechazo));
        
        serie_ActHjorth = removerows(serie_ActHjorth_sin_rechazo','ind',find(eliminar_latidos_aux));
        serie_Hilbert = removerows(serie_Hilbert_sin_rechazo','ind',find(eliminar_latidos_aux));
        t = removerows(t_sin_rechazo','ind',find(eliminar_latidos_aux));
        
        if d==1
            t_sin_rechazo1 = t_sin_rechazo';
            serie_ActHjorth_sin_rechazo1 = serie_ActHjorth_sin_rechazo';
            [serie_ActHjorth1,t1] = filtrado_y_normalizado(serie_ActHjorth,t,100);
            serie_Hilbert_sin_rechazo1 = serie_Hilbert_sin_rechazo';
            [serie_Hilbert1,~] = filtrado_y_normalizado(abs(serie_Hilbert),t,100);
            eliminar_latidos1 = eliminar_latidos_aux;
            
        else
            t_sin_rechazo2 = t_sin_rechazo';
            serie_ActHjorth_sin_rechazo2 = serie_ActHjorth_sin_rechazo';
            [serie_ActHjorth2,t2] = filtrado_y_normalizado(serie_ActHjorth,t,100);
            serie_Hilbert_sin_rechazo2 = serie_Hilbert_sin_rechazo';
            [serie_Hilbert2,~] = filtrado_y_normalizado(abs(serie_Hilbert),t,100);
            eliminar_latidos2 = eliminar_latidos_aux;
        end
       
        clearvars serie_ActHjorth_sin_rechazo t_sin_rechazo serie_Hilbert_sin_rechazo
    end
    
%% unifico las dos series en una, promediando punto a punto
    if length(t1) > length(t2)
        serie_ActHjorth_Unif = sqrt(0.5*((serie_ActHjorth1(1:length(t2))).^2+(serie_ActHjorth2(1:length(t2))).^2));
        serie_Hilbert_Unif = sqrt(0.5*((serie_Hilbert1(1:length(t2))).^2+(serie_Hilbert2(1:length(t2))).^2));
        t_Unif = t2;                                                        
    else
        serie_ActHjorth_Unif = sqrt(0.5*((serie_ActHjorth1(1:length(t1))).^2+(serie_ActHjorth2(1:length(t1))).^2));
        serie_Hilbert_Unif = sqrt(0.5*((serie_Hilbert1(1:length(t1))).^2+(serie_Hilbert2(1:length(t1))).^2));
        t_Unif = t1;
    end
    
    save([direcGuardado '\' reg],'eliminar_latidos','serie_ActHjorth1','serie_ActHjorth2','serie_ActHjorth_sin_rechazo1','serie_ActHjorth_sin_rechazo2','t1','t2','t_sin_rechazo1','t_sin_rechazo2','serie_Hilbert1','serie_Hilbert2','serie_Hilbert_sin_rechazo1','serie_Hilbert_sin_rechazo2','serie_ActHjorth_Unif','t_Unif','serie_Hilbert_Unif','eliminar_latidos1','eliminar_latidos2','ruido')
end