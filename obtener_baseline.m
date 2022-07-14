function baseline = obtener_baseline(ECG,Poff,QRSon,inicio,fin,flag_plot)
% Toma la señal, le calcula la línea de base (usando cubic spline) y  
% devuelve la linea de base

%% entrada:
%   ECG: señal de nx1 
%   Poff: fin de la onda P para cada latido
%   QRSon: comienzo del complejo QRS para cada latido
%   inicio: comienzo para calcular la base de línea (ej: 1)
%   fin: fin para calcular la base d elínea (ej: 1800000)
%   flag_plot: 1 o 0 si se quiere graficar los resultados obtenidos

%% salida
%   baseline: línea de base

%%
%% Busco el primer y último valor de las marcas de P y R que estén dentro del intervalo [inicio;fin] que estoy analizando

inicio_PR = find(Poff>=inicio);
inicio_PR = inicio_PR(1);
fin_PR = find(QRSon>=fin);
if isempty(fin_PR)
    fin_PR = length(Poff);
else
    fin_PR = fin_PR(1); 
end
    
if(flag_plot)
    figure;
    plot(ECG(inicio:fin)); % señal original
    hold on;
end

%% Busco los valores medios del segmento PR ("knots")

PRmed = zeros(fin_PR-inicio_PR+1+2,1);
PRmed(1) = inicio;
PRmed(length(PRmed(:,1))) = fin;
PRmed_amp = zeros(fin_PR-inicio_PR+1+2,1);
PRmed_amp(1) = ECG(inicio);
PRmed_amp(length(PRmed(:,1))) = ECG(fin);


desv = 0;
aux=[];
for i = inicio_PR:fin_PR
    if ~isnan(Poff(i))
        if ~isnan(QRSon(i))
             PRmed(i-inicio_PR+2,1) = round( (Poff(i)+QRSon(i)) / 2.0 );
             desv = PRmed(i-inicio_PR+2,1) - Poff(i);
        else
             PRmed(i-inicio_PR+2,1) = Poff(i) + desv;
        end
    elseif ~isnan(QRSon(i))
         PRmed(i-inicio_PR+2,1) = QRSon(i) - desv;
    else
        PRmed(i-inicio_PR+2,1) = PRmed(i-1-inicio_PR+2,1);
        aux=[aux i-inicio_PR+2];
    end
    
    PRmed_amp(i-inicio_PR+2,1) = ECG(PRmed(i-inicio_PR+2,1));
    
    if(flag_plot)
         scatter(PRmed(i-inicio_PR+2,1),PRmed_amp(i-inicio_PR+2,1),3,'MarkerEdgeColor','red');
    end
end

PRmed = removerows(PRmed,aux);
PRmed_amp = removerows(PRmed_amp,aux);


%% Calculo la línea de base haciendo un cubic spline con los PRmed

baseline = interp1(PRmed(:,1),PRmed_amp(:,1),inicio:fin,'spline');

%% Grafico la línea de base
if(flag_plot)
    plot(baseline,'color','red');
    hold off;
end