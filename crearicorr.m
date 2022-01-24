function nodo = crearicorr(FormaPolino,Th,PuntosGeo,e,l)
% Función que permite crear un fichero que guarda la intensidad que
% maximiza la diferencia de temperaturas entre la temperatura en el foco
% caliente y el foco frío para diferentes formas de termoelemento.

% Regresión de ImxAp con parámetros geométricos
if Th==50+273.15 && e==1e-3 && FormaPolino==3 && l>=2.5e-3 && l<= 7e-3 && max(PuntosGeo)<=7.5e-4 && min(PuntosGeo)>=3.75e-4
    load RegresionVar
    ImxAp = interpn(Param1,Param2,Param3,Param4,intensi,l,PuntosGeo(1),PuntosGeo(2),PuntosGeo(3),'spline');
else
    a0 = 1.988e-4;
    a1 = 3.53e-7;
    a2 = 7.521e-10;
    g0 = 1.096e5;
    g1 = -5.59e2;
    g2 = 2.498;
    k0 = 1.663;
    k1 = -3.58e-3;
    k2 = 3.195e-5;
    alpha = 1.954e-4;
    gamma = 1.154e5;
    kappa = 1.702;
    
    Tcest0 = (-kappa + sqrt(kappa*(2*Th*alpha^2*gamma + kappa))) / (alpha^2*gamma);
    
    error = realmax;
    tol   = 1e-12;
    while error>tol
        Tmed = (Tcest0 + Th)/2 - 273.15;
        alpha = a0 + a1 * Tmed + a2 * Tmed^2;
        gamma = g0 + g1 * Tmed + g2 * Tmed^2;
        kappa = k0 + k1 * Tmed + k2 * Tmed^2;
        Tcest1 = (-kappa + sqrt(kappa*(2*Th*alpha^2*gamma + kappa))) / (alpha^2*gamma);
        error  = abs(Tcest1-Tcest0)/Tcest0;
        Tcest0 = Tcest1;
    end
    
    if FormaPolino == 2
        d1 = 2*e/l*(PuntosGeo(2)-PuntosGeo(1));
        d2 = 2*e*PuntosGeo(1);
        if PuntosGeo(1) == PuntosGeo(2)
            ImxAp = d2*Tcest0*alpha*gamma/l;
        else
            ImxAp = (-d1*Tcest0*alpha*gamma)/(log(d2*gamma*kappa)-log((d2+d1*l)*gamma*kappa));
        end
    elseif FormaPolino == 3
        d1 = 4*e/l^2*(PuntosGeo(1)-2*PuntosGeo(2)+PuntosGeo(3));
        d2 = -2*e/l*(3*PuntosGeo(1)-4*PuntosGeo(2)+PuntosGeo(3));
        d3 = 2*e*PuntosGeo(1);
        raiz = sqrt(-d2^2+4*d1*d3);
        if PuntosGeo(1) == PuntosGeo(2) && PuntosGeo(2) == PuntosGeo(3)
            ImxAp = d3*Tcest0*alpha*gamma/l;
        elseif PuntosGeo(2) == ((PuntosGeo(1)+PuntosGeo(3))/2)
            ImxAp = (-d2*Tcest0*alpha*gamma)/(log(d3*gamma*kappa)-log((d3+d2*l)*gamma*kappa));
        else
            ImxAp = (-raiz*Tcest0*alpha*gamma)/(2*(atan(d2/raiz)-atan((d2+2*d1*l)/raiz)));
        end
    end
end

% Modificación del número de elementos con la altura.
nT = ceil(chop(l/3e-4,6))*2;
nodo = 573+30*(nT/2-1);

if FormaPolino == 2
    fid=fopen('icorr','w');
    fprintf(fid,'PARAmeter \n');
    fprintf(fid,'jz = %e \n',ImxAp);
    fprintf(fid,'L1 = %e \n',PuntosGeo(1));
    fprintf(fid,'L2 = %e \n',PuntosGeo(2));
    fprintf(fid,'nf = %i \n',nodo);
    fprintf(fid,'nT = %i \n\n',nT);
    fclose(fid);
elseif FormaPolino == 3
    fid=fopen('icorr','w');
    fprintf(fid,'PARAmeter \n');
    fprintf(fid,'jz = %e \n',ImxAp);
    fprintf(fid,'L1 = %e \n',PuntosGeo(1));
    fprintf(fid,'L2 = %e \n',PuntosGeo(2));
    fprintf(fid,'L3 = %e \n',PuntosGeo(3));
    fprintf(fid,'nf = %i \n',nodo);
    fprintf(fid,'nT = %i \n\n',nT);
    fclose(fid);
end
