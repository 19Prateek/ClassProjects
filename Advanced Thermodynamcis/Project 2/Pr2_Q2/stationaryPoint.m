function [y,Xi,iter] = stationaryPoint(Xi,Zfixed,b,a,B,aij,Kbinary,Tc,Rgas,Temp,P,ncomps,Z,phiZfixed)
lnphi=Tc;f=Tc;phi=Tc;
y=[0,0];
maxIter=500;
iter=0;
tol=1e-10;
criteria=Tc;

while norm(criteria) > tol && (iter < maxIter) 
    y=Xi/sum(Xi);  
    %Mixture properties : bmix
    bmix=0;
    for i=1:ncomps
        bmix=y(i)*b(i) + bmix;
    end
    %Mixture properties : amix
    amix=0;
    for i=1:ncomps
        for j=1:ncomps
            aij(i,j)=sqrt(a(i)*a(j))*(1-Kbinary(i,j));
            amix=amix+y(i)*y(j)*aij(i,j);
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
            Zmix=Zlmix;
        elseif delGRTvmix<delGRTlmix
            Zmix=Zvmix;
        end
    else
        Zmix=Zmixture;
    end
    %Fugacity Equation for each species in mixture
    for i=1:ncomps
        sumterm=0;
        for j=1:ncomps
            sumterm= sumterm + y(j)*Aij(i,j);
        end
        lnphi(i)=((B(i)/Bmix)*(Zmix-1))-log(Zmix-Bmix)-(((2*sumterm)/Amix)-(B(i)/Bmix))*(Amix/(2*sqrt(2)*Bmix))*log((Zmix+(1+sqrt(2))*Bmix)/(Zmix+(1-sqrt(2))*Bmix));
        phi(i)=exp(lnphi(i));
        f(i)=phi(i)*y(i)*P/6894.76;
    end     
    criteria =abs(log(Xi.*phi)-log(Zfixed.*phiZfixed));
    Xi=(Zfixed.*phiZfixed)./(phi);
    iter = iter +1;
end
