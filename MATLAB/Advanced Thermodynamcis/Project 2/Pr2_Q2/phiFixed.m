function [ phiZfixed,Zfixed] = phiFixed(Zfixed,b,a,B,Kbinary,aij,P,Rgas,Temp,Tc,ncomps,Z)

%******************Stability Analysis : fixed Composition calculation******
phiZfixed=Tc;lnphifixed=Tc;
yfixed=Zfixed;
%Mixture properties : bmix
bmix=0;
for i=1:ncomps
    bmix=yfixed(i)*b(i) + bmix;
end
%Mixture properties : amix
amix=0;
for i=1:ncomps
    for j=1:ncomps
        aij(i,j)=sqrt(a(i)*a(j))*(1-Kbinary(i,j));
        amix=amix+yfixed(i)*yfixed(j)*aij(i,j);
    end
end
Aij=(aij*P)/((Rgas*Temp)^2);
Amix=(amix*P)/((Rgas*Temp)^2);
Bmix=(bmix*P)/(Rgas*Temp);

%Cordano Cubic Roots for Peng-Robinson
[~,coeff]=cordano(Amix,Bmix,Z);
%Cubic Roots for Compressibility factor
[Zmixture,nbr_roots]=cubic_roots(coeff);
if nbr_roots>1
    Zlmix=min(Zmixture);
    Zvmix=max(Zmixture);
    delGRTlmix=(Zlmix-1)-log(Zlmix-Bmix)-(Amix/(2*sqrt(2)*Bmix))*log((Zlmix+(1+sqrt(2))*Bmix)/(Zlmix+(1-sqrt(2))*Bmix)) ;
    delGRTvmix=(Zvmix-1)-log(Zvmix-Bmix)-(Amix/(2*sqrt(2)*Bmix))*log((Zvmix+(1+sqrt(2))*Bmix)/(Zvmix+(1-sqrt(2))*Bmix)) ;
    if delGRTlmix<delGRTvmix
        Zmix=min(Zmixture);
    else
        Zmix=max(Zmixture);
    end
else
    Zmix=Zmixture;
end

%Fugacity Equation for each species in mixture
for i=1:ncomps
    sumterm=0;
    for j=1:ncomps
        sumterm= sumterm + yfixed(j)*Aij(i,j);
    end
    lnphifixed(i)=((B(i)/Bmix)*(Zmix-1))-log(Zmix-Bmix)-(((2*sumterm)/Amix)-(B(i)/Bmix))*(Amix/(2*sqrt(2)*Bmix))*log((Zmix+(1+sqrt(2))*Bmix)/(Zmix+(1-sqrt(2))*Bmix));
    phiZfixed(i)=exp(lnphifixed(i));
end

end

