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

% Registros que voy a descomponer (desde primer_registro hasta ultimo_registro)
primer_registro = 'e0103'; 
ultimo_registro = 'e1304';
in = strfind(record,primer_registro);
fi = strfind(record,ultimo_registro);

carpeta_guardado = '<file_path>'; % <file_path> debe ser cambiado por la dirección donde se van a guardar los resultados del EMD

% Recorro cada registro
for r = in : 5 : fi
    reg = record(r:r+4);
    disp(reg);
    
    load('<file_path>') % <file_path> debe ser cambiado por la dirección donde se encuentran los registros de la base de datos en .mat (variable: registro)
    
    mkdir([carpeta_guardado '\' reg]);
    
    for d=1:2
        ecg = registro.ECG(:,d+1);
        T = registro.ECG(:,1);
        fs = registro.FS;
            
        [imf,~] = emd(ecg);
        
        save('<file_path>') % <file_path> debe ser cambiado por la dirección donde se guardaran los resultados
    end
end
