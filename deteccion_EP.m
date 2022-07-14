function [EP_flag,inicio_EP,fin_EP,baseline,umbral] = deteccion_EP(serie,t_muestras,beta,nu)
% A partir del beta y nu indicado detecta episodios isquémicos de la señal
% utilizando la serie, devuelve un vector de flags en donde se indica
% con 1 si hay un episodio isquémico en ese latido o 0 si no hay. También
% devuelve dos vectores con los momentos en el que comienzan y terminan los
% episodios

%% entrada:
%   serie: serie obtenida a partir de los parámetros
%   t_muestras: vector de tiempo con número de muestras
%   beta: variable que determina cuanto se ajusta la línea de base a la
%         serie, para luego detectar alteraciones (episodios isquémicos)
%   nu: variable que determine en cuanto tiene que superar la serie a
%       la línea de base para considerarlo un episodio isquémico

%% salida
%   EP_flag: vector del largo de la cantidad de latidos, en dónde se indica
%            con 1 o 0 según si hay un episodio isquémico o no, respectivamente
%   inicio_EP: vector con los comienzos de los episodios isquémicos
%   fin_EP: vector con el momento en el que terminan los episodios isquémicos
%   baseline: línea de base de la serie, generada a partir del beta
%   umbral: umbral del detector

%%
%%
EP_flag = zeros(length(serie),1);

baseline = zeros(length(serie),1);
baseline(1) = serie(1);
cant_EP = 0;
aux = 0;
umbral = zeros(length(serie),1);
umbral(1) = baseline(1) + nu;

% Punto a punto de la serie se va calculando la línea de base, sólo se
% calcula cuando no se considera que es parte de un episodio isquémico (EP)
% Se considera que es un EP cuando la serie supera al umbral, que es la
% línea de base + nu

for n=2:length(serie)
    
    umbral(n) = baseline(n-1) + nu;
    if serie(n) > umbral(n)
        EP_flag(n) = 1;
        baseline(n) = baseline(n-1);

        if EP_flag(n-1) == 0 
            cant_EP = cant_EP + 1;
            inicio_EP(cant_EP) = t_muestras(n);
            inicio_EP_n(cant_EP) = n;
        end
    else
        if EP_flag(n-1) == 0
            baseline(n) = baseline(n-1) + beta * (serie(n)-baseline(n-1));
            aux = aux + 1;
        else
            fin_EP(cant_EP) = t_muestras(n);
            fin_EP_n(cant_EP) = n;
            duracion_EP = fin_EP(cant_EP)-inicio_EP(cant_EP);
            cant_latidos_EP = fin_EP_n(cant_EP)-inicio_EP_n(cant_EP);
            
            if duracion_EP<45.0*250 || cant_latidos_EP<4
                EP_flag(inicio_EP_n(cant_EP):fin_EP_n(cant_EP)) = 0;
                inicio_EP = inicio_EP(1:cant_EP-1);
                fin_EP = fin_EP(1:cant_EP-1);
                inicio_EP_n = inicio_EP_n(1:cant_EP-1);
                fin_EP_n = fin_EP_n(1:cant_EP-1);
                cant_EP = cant_EP - 1;
            end
            baseline(n) = baseline(n-1);
        end
    end
end

if exist('inicio_EP','var') && exist('fin_EP','var')
    if length(inicio_EP) > length(fin_EP)
        fin_EP(cant_EP) = t_muestras(n);
        fin_EP_n(cant_EP) = n; 
        duracion_EP = fin_EP(cant_EP)-inicio_EP(cant_EP);
        if duracion_EP<45.0*250
            EP_flag(inicio_EP_n(cant_EP):fin_EP_n(cant_EP)) = 0;
            inicio_EP = inicio_EP(1:cant_EP-1);
            fin_EP = fin_EP(1:cant_EP-1);
            inicio_EP_n = inicio_EP_n(1:cant_EP-1);
            fin_EP_n = fin_EP_n(1:cant_EP-1);
            cant_EP = cant_EP - 1;
        end
    end
end

    
if exist('inicio_EP','var')
    if ~exist('fin_EP','var')
        fin_EP = length(serie);
    end
    
    a=0;
    for i=1:length(inicio_EP-1)
        if i-a >= length(inicio_EP)
            break;
        elseif inicio_EP(i+1-a)-fin_EP(i-a)<120*250
            EP_flag(fin_EP_n(i-a):inicio_EP_n(i+1-a)) = 1;
            inicio_EP = removerows(inicio_EP',i+1-a);
            inicio_EP = inicio_EP';
            inicio_EP_n = removerows(inicio_EP_n',i+1-a);
            inicio_EP_n = inicio_EP_n';
            fin_EP(i-a) = fin_EP(i+1-a);
            fin_EP_n(i-a) = fin_EP_n(i+1-a);
            fin_EP = removerows(fin_EP',i+1-a);
            fin_EP = fin_EP';
            fin_EP_n = removerows(fin_EP_n',i+1-a);
            fin_EP_n = fin_EP_n';
            cant_EP = cant_EP - 1;
        end
    end
end

if cant_EP == 0
    inicio_EP = 0;
    fin_EP = 0;
end

