function [iter3,K,y,x,V,L] = TwoPFlash(Ki,Pc,Tc,w,P,Kbinary,Temp,Rgas,ncomps,aij,Zfixed,Z)

K=Ki;
%disp (['K_Wilson =    ' num2str(K) ]);
maxnorm=Tc;
phiZV=Tc; phiZL=Tc;
lnphiV=Tc;lnphiL=Tc;
iter3=0;
Zi=Zfixed;
maxIter=2000;
tol=1e-10;

%RR Equation
while (norm(maxnorm)>tol) && (iter3 < maxIter)
    V=-3:0.0001:1.5;
    fV=0;
    fVMAT=V;
    for j=1:size(V,2)
        for i=1:ncomps
            fV = fV + ((1-K(i))*Zi(i))/ (1-((1-K(i))*V(j)));
        end
        fVMAT(j)=fV;
        fV=0;
    end
    Kmin=min(K);
    Kmax=max(K);
    PolesCorrect=[1/(1-Kmax) 1/(1-Kmin)];
    Poles=1./(1-K);
    
    %Newton Method to solve for V
    delV=0.000001;
    Vnew=1;
    Vold=0;
    fVMAT=Zi;
    fVdelMAT=Zi;
    while abs(Vnew-Vold) > tol
        Vold=Vnew;
        for i=1:ncomps
            fVMAT(i) =  ((1-K(i))*Zi(i))/ (1-((1-K(i))*Vold));
            fV= sum(fVMAT);
            fVdelMAT(i) =  ((1-K(i))*Zi(i))/ (1-((1-K(i))*(Vold+delV)));
            fVdel =sum(fVdelMAT);
        end
        derfVold = (fVdel-fV)/delV;
        Vnew= Vold - (fV/derfVold);
    end
    %Phase Fractions
    V=Vnew;
    L=1-V;
    %Component Fractions
    y=Zi;x=Zi;
    for i=1:ncomps
        x(i)=Zi(i)/(L+K(i)*(1-L));
        y(i)=K(i)*x(i);
    end
%        
%     disp (['Vapour Fractions =    ' num2str(y) ]);
%     disp (['Liquid Fractions =    ' num2str(x) ]);
%     disp (['V at Poles =    ' num2str(Poles) ]);
%     disp (['Correct Solution Exists between =    ' num2str(PolesCorrect) ]);
%     disp (['V =  ' num2str(V)]);
%     disp (['L =  ' num2str(L)]);
    
    
   %******************************VAPOUR PHASE*****************************
    
    %Calculation of cubic EOS parameters
    [a,b,A,B]=PR_para(Pc,Tc,w,P,Temp,Rgas,ncomps);
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
    [ZV,nbr_roots]=cubic_roots(coeff);
    if nbr_roots>1
        Zl=min(ZV);
        Zv=max(ZV);
        delGRTl=(Zl-1)-log(Zl-B(i))-(A(i)/(2*sqrt(2)*B(i)))*log((Zl+(1+sqrt(2))*B(i))/(Zl+(1-sqrt(2))*B(i))) ;
        delGRTv=(Zv-1)-log(Zv-B(i))-(A(i)/(2*sqrt(2)*B(i)))*log((Zv+(1+sqrt(2))*B(i))/(Zv+(1-sqrt(2))*B(i))) ;
        if delGRTl<delGRTv 
            ZmixV=Zl;
        elseif delGRTv<delGRTl 
            ZmixV=Zv;
        end
    else
        ZmixV=ZV;
    end
    
    %Fugacity Equation for each species in mixture
    for i=1:ncomps
        sumterm=0;
        for j=1:ncomps
            sumterm= sumterm + y(j)*Aij(i,j);
        end
        lnphiV(i)=((B(i)/Bmix)*(ZmixV-1))-log(ZmixV-Bmix)-(((2*sumterm)/Amix)-(B(i)/Bmix))*(Amix/(2*sqrt(2)*Bmix))*log((ZmixV+(1+sqrt(2))*Bmix)/(ZmixV+(1-sqrt(2))*Bmix));
        phiZV(i)=exp(lnphiV(i));
    end
    
    
    %******************************LIQUID PHASE*****************************
    %Mixture properties : bmix
    bmixL=0;
    for i=1:ncomps
        bmixL=x(i)*b(i) + bmixL;
    end   
    %Mixture properties : amix
    amixL=0;
    for i=1:ncomps
        for j=1:ncomps
            aij(i,j)=sqrt(a(i)*a(j))*(1-Kbinary(i,j));
            amixL=amixL+x(i)*x(j)*aij(i,j);
        end
    end
    Aij=(aij*P)/((Rgas*Temp)^2);
    AmixL=(amixL*P)/((Rgas*Temp)^2);
    BmixL=(bmixL*P)/(Rgas*Temp);
    
    %Cordano Cubic Roots for Peng-Robinson
    [~,coeff]=cordano(AmixL,BmixL,Z);
    %Cubic Roots for Compressibility factor
    [ZL,nbr_roots]=cubic_roots(coeff);
    if nbr_roots>1
        Zl=min(ZL);
        Zv=max(ZL);
        delGRTl=(Zl-1)-log(Zl-B(i))-(A(i)/(2*sqrt(2)*B(i)))*log((Zl+(1+sqrt(2))*B(i))/(Zl+(1-sqrt(2))*B(i))) ;
        delGRTv=(Zv-1)-log(Zv-B(i))-(A(i)/(2*sqrt(2)*B(i)))*log((Zv+(1+sqrt(2))*B(i))/(Zv+(1-sqrt(2))*B(i))) ;
        if delGRTl<delGRTv 
            ZmixL=Zl;
        elseif delGRTv<delGRTl 
            ZmixL=Zv;
        end
    else
        ZmixL=ZL;
    end
    
    %Fugacity Equation for each species in mixture
    for i=1:ncomps
        sumterm=0;
        for j=1:ncomps
            sumterm= sumterm + x(j)*Aij(i,j);
        end
        lnphiL(i)=((B(i)/BmixL)*(ZmixL-1))-log(ZmixL-BmixL)-(((2*sumterm)/AmixL)-(B(i)/BmixL))*(AmixL/(2*sqrt(2)*BmixL))*log((ZmixL+(1+sqrt(2))*BmixL)/(ZmixL+(1-sqrt(2))*BmixL));
        phiZL(i)=exp(lnphiL(i));
    end 
    maxnorm=abs(log(x.*phiZL)-log(y.*phiZV));
    K=exp(log(phiZL)-log(phiZV));
    iter3=iter3+1;
end
%     disp (['Vapour Fractions =    ' num2str(y) ]);
%     disp (['Liquid Fractions =    ' num2str(x) ]);
%     disp (['V at Poles =    ' num2str(Poles) ]);
%     disp (['Correct Solution Exists between =    ' num2str(PolesCorrect) ]);
%     disp (['V =  ' num2str(V)]);
%     disp (['L =  ' num2str(L)]);
end

