% ------------ SCRIPT PRINCIPAL PARA 'SIMULATED ANNEALING'  ------------ %
% ------------ APLICADO A OPTIMIZACION DE TERMOELECTRICOS   ------------ %

% Cargar Media y Variación típica del caso.

load ParamNorm

% Definicion de intervalos maximos

IntAncho    = [2,10];             % Intervalo de valores de ancho de pulso
IntLongitud = [2.5,7]*1e-3;       % Intervalo de valores de longitud de termoelemento
IntAnchoGeo = [0.75,1.5]*1e-3/2;  % Intervalo de valores de ancho geométrico
NumPuntPuls = 2;                  % Número de puntos de pulso
sigma_adm   = 60e6;               % Tensión máxima admisible
FormaPolino = 3;				  % 2 es lineal, 3 es cuadrático
Th          = 50+273.15;
e           = 1e-3;

% Definicion de parametros de la función de Boltzmann

t0_a = 0.12;                              % Temperatura inicial
lambda_a = 0.85;                            % Parametro de reduccion de temperatura (0 < lambda < 1)
Kmax = 100*(1+FormaPolino+1);               % Unas 100 veces el numero de parametros
Kaccept = 10*(1+FormaPolino+1);             % Unas 10 veces el numero de parametros
Stopping = 3;                               % Numero para parar el algoritmo
Almacen = 6*Kmax;                           % Inicialización de ciertas variables
load alfa
alfa   = alfa'/norm(alfa);

% Parametros iniciales a introducir

AnchoPul = mean(IntAncho);                  % Pulso de duración media
Longitud = mean(IntLongitud);               % La longitud del termoelemento es la media del intervalo.
GananciaIP = 2.5*ones(1,NumPuntPuls);       % Pulso cuadrado, de ganancia 2,5.
PuntosGeo = ones(1,FormaPolino)*1e-3/2;     % Forma constante.

% Preparar ipulse

q        = 500;                             % Máximos instantes de tiempo calculados en cada caso.
dtc      = 0.1;                             % Paso temporal en la parte constante de la resolución del caso
AnchoPul = AnchoPul + dtc*(NumPuntPuls-1) - mod(AnchoPul,dtc*(NumPuntPuls-1)); % Hace que el Ancho de Pulso sea divisible por dtc*NumPuntPuls
HisParam = zeros(1+FormaPolino+1+7,Almacen);% Variable que acumula los pulsos considerados y la determinación tomada.
tiem     = zeros(q,Almacen);                % Acumula todos los pasos de tiempo de cada uno de los casos corridos
temp     = zeros(q,Almacen);                % Acumula la historia de la temperatura de cada uno de los casos corridos
volt     = zeros(q,Almacen);                % Acumula la historia del voltaje de cada uno de los casos corridos
tens     = zeros(1,Almacen);                % Acumula máxima tensión de Tresca de cada uno de los casos corridos

% Parametros de bucle

NombreFichero = 'if';                       % El nombre del fichero de FEAP que tiene que correr.
acabar   = 0;                               % Variable que es 1 cuando ha llegado al final de la optimización
iter     = 0;                               % La iteración
PARAM    = [AnchoPul,PuntosGeo,Longitud];            % Parámetros iteracion previa
PARAMN   = [AnchoPul,PuntosGeo,Longitud];            % Parámetros nueva Iteracion
m        = length(PARAM);                   % 
O1       = realmax;                         % O1: Función de coste. Primer valor infinito.
t        = t0_a;                            % temperatura de annealing. Inicializada con t0_a
iteracep = 0;                               % Número de iteraciones aceptadas
iterrech = 0;                               % Número de iteraciones rechazadas
numKmax  = 0;                               % Número de veces que se alcanza Kmax en iterrech.
x2       = 0;                               % Variable que acumula la historia temporal de las variables temperatura y voltaje.

%Inicializar rand; con un valor de entrada

%% Bucle Principal ---- NO TOCAR NADA. PARAMETROS ARRIBA

while ~acabar
    % CALCULO DEL CASO
    iter = iter + 1;
    HisParam(1:m,iter) = PARAMN;
    crearicorr(FormaPolino,Th,PARAMN(1+(1:FormaPolino)),e,PARAMN(m));
    inputsfeap(NombreFichero,PARAMN(m));
    crearipulse(GananciaIP,PARAMN(1),dtc);
    vx        = unix(['feap -i',NombreFichero]);
    fprintf('Acabado el caso %i \n',iter)
    x         = load(['P',NombreFichero(2:length(NombreFichero)),'a.dis']);
    b         = length(x);
    x(b+1:q,:) = nan;
    tiem(:,iter) = x(:,1);
    temp(:,iter) = x(:,2);
    volt(:,iter) = x(:,3)-x(:,4);
    
    % HALLAR MÍNIMO DE TEMPERATURA. CON PRECISIÓN.
    [~,pos1]    = min(temp(1:b,iter));
    vectiem     = linspace(tiem(pos1-1,iter),tiem(pos1+1,iter),10000);
    [Tmin,pos2] = min(interp1(tiem(1:b,iter),temp(1:b,iter),vectiem,'spline'));
    Tss1        = temp(1,iter);
    DeltaTp     = Tss1 - Tmin;
    [~,pos3]    = max(temp(1:b,iter));
    vectiem2    = linspace(tiem(pos3-1,iter),tiem(pos3+1,iter),10000);
    Tmax        = max(interp1(tiem(1:b,iter),temp(1:b,iter),vectiem2,'spline'));
    DeltaTpp    = Tmax - Tss1;
    vecHolding  = temp(1:b,iter)+abs(Tss1)+DeltaTp*0.8 < 0;
    posicionh1  = find(diff(vecHolding)>0) + 1;
    posicionh2  = find(diff(vecHolding)<0);
    longposih   = length(posicionh1);
    th1        = zeros(1,longposih);
    th2        = zeros(1,longposih);
    for i=1:longposih
        vectiem3    = linspace(tiem(posicionh1(i)-1,iter),tiem(posicionh1(i)+1,iter),10000);
        vectiem4    = linspace(tiem(posicionh2(i)-1,iter),tiem(posicionh2(i)+1,iter),10000);
        [~,pos4]    = min(interp1(tiem(1:b,iter),abs(temp(1:b,iter)+abs(Tss1)+DeltaTp*0.8),vectiem3,'spline'));
        [~,pos5]    = min(interp1(tiem(1:b,iter),abs(temp(1:b,iter)+abs(Tss1)+DeltaTp*0.8),vectiem4,'spline'));
        th1(i)      = vectiem3(pos4);
        th2(i)      = vectiem4(pos5);
    end
    
    
    f1 = 1/DeltaTp;                 % Temperatura mínima
    f2 = vectiem(pos2)-dtc;         % Tiempo tmin
    f3 = DeltaTpp;                  % Temperatura maxima
    f4 = 1/(sum(th2-th1));          % Holding time
    
    O2 = alfa*(([f1,f2,f3,f4]-[f1media,f2media,f3media,f4media])./[f1destip,f2destip,f3destip,f4destip])';
    HisParam(m+2+(1:5),iter) = [f1,f2,f3,f4,O2];
    
    
    tresca = load('tresca');
    tens(iter) = tresca(1);
    if tresca(1) < sigma_adm
        
        if O2 < O1
            fprintf('El caso mejora la función de coste. O1 - O2 = %g\n',O1-O2)
            PARAM = PARAMN;
            iteracep = iteracep + 1;
            O1 = O2;
            x2 = x;
            HisParam(m+1,iter) = 1;
            
        else
            if O2 == O1 && j<=NumPuntPuls % Si la ganancia de un punto no influye, se toma la menor.
                if PARAMN(j) > PARAM(j)
                    PARAMN(j) = PARAM(j);
                end
            end
            PBoltzmann = exp((O1-O2)/t);
            l = rand(1);
            if l <= PBoltzmann
                fprintf('El caso no mejora la función de coste, pero se acepta. O1 - O2 = %g\n',O1-O2)
                PARAM = PARAMN;
                iteracep = iteracep + 1;
                O1 = O2;
                x2 = x;
                HisParam(m+1,iter) = 2;
            else
                fprintf('El caso no mejora la función de coste y no se acepta. O1 - O2 = %g\n',O1-O2)
                PARAMN = PARAM;
                iterrech = iterrech + 1;
                HisParam(m+1,iter) = 3;
            end
            
        end
        
    else
        fprintf('Las tensiones máximas superan las admisibles. El caso se rechaza. O1 - O2 = %g\n',O1-O2)
        PARAMN = PARAM;
        iterrech = iterrech + 1;
        HisParam(m+1,iter) = 4;
    end
    
    
    if iteracep >= Kaccept
        t = lambda_a * t;
        iteracep = 0;
    end
    
    if iterrech >= Kmax
        t = lambda_a * t;
        numKmax = numKmax + 1;
        iterrech = 0;
    end
    
    if numKmax == Stopping
        acabar = 1;
        save ResultadosFinales
    end
    
    % Siguiente iteracion, preparacion de parametros
    
    j = ceil(rand(1)*m);        % Parametro a cambiar
    HisParam(m+2,iter) = j;
    if ~acabar
        if j == m
            PARAMN(j) = IntLongitud(1)+rand(1)*(IntLongitud(2)-IntLongitud(1));
        elseif j == 1
            PARAMN(j) = IntAncho(1)+rand(1)*(IntAncho(2)-IntAncho(1));
            PARAMN(j) = PARAMN(j) + dtc*(NumPuntPuls-1) - mod(PARAMN(j),dtc*(NumPuntPuls-1));
            while PARAMN(j) == PARAM(j)
                PARAMN(j) = IntAncho(1)+rand(1)*(IntAncho(2)-IntAncho(1));
                PARAMN(j) = PARAMN(j) + dtc*(NumPuntPuls-1) - mod(PARAMN(j),dtc*(NumPuntPuls-1));
            end
        else
            PARAMN(j) = IntAnchoGeo(1)+rand(1)*(IntAnchoGeo(2)-IntAnchoGeo(1));
        end
    end
    
    if ~mod(iter,5)
        save('resultados','PARAM','PARAMN','iteracep','iterrech','numKmax','x2','HisParam','O1','f1','f2','f3','f4','t','iter','tens','tiem','temp','volt')
    end
    
    
end
