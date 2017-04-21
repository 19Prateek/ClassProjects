%clear;
% close all;
clc;
%Properties

Z=[0,0,0];
% Component 1:CO2, Component 2:nC4, Component 3:nC7, Component 4:nC
Temp = 621.8 ;                  %K
P=35.49e5;                     %Pa
Tc= [507.60 871.16 ];       %K
Pc= [30.25e5 5.53e5 ];     %Pa
w= [0.3010 1.4678];
Rgas= 8.3144598;  %Pa m^3 K^-1 mol^-1
ncomps=2;
Zi=[0 0];
K =[0 0;
    0 0];
aij=K;
lnphi=Tc;
f=Tc;
phi=Tc;
phiPure=Tc;
phiZfixed=Tc;
lnphifixed=Tc;
%**********************************
%Peng Robinson EOS parameters
[a,b,A,B,k,alpha]=PR_para(Pc,Tc,w,P,Temp,Rgas,ncomps);

xC1=0:0.001:0.999;
delGmixRTMAT=xC1;
delGmixRT1MAT=xC1;
DRMAT=xC1;
phiMAT=zeros(size(xC1,2),ncomps);
phiPureMAT=phiMAT;

%*******************************fixed Composition calculation**************
Zfixed= [1-0.0502 0.0502];
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
        aij(i,j)=sqrt(a(i)*a(j))*(1-K(i,j));
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
sumterm=0;
for i=1:ncomps
    sumterm=0;
    for j=1:ncomps
        sumterm= sumterm + yfixed(j)*Aij(i,j);
    end
    lnphifixed(i)=((B(i)/Bmix)*(Zmix-1))-log(Zmix-Bmix)-(((2*sumterm)/Amix)-(B(i)/Bmix))*(Amix/(2*sqrt(2)*Bmix))*log((Zmix+(1+sqrt(2))*Bmix)/(Zmix+(1-sqrt(2))*Bmix));
    phiZfixed(i)=exp(lnphifixed(i));
end

%**************************************************************************
for counter=1:size(xC1,2)
    Zi=[xC1(counter) 1-xC1(counter)];
    
    %Mixture properties : bmix
    y=Zi;
    bmix=0;
    for i=1:ncomps
        bmix=y(i)*b(i) + bmix;
    end
    
    %Mixture properties : amix
    amix=0;
    for i=1:ncomps
        for j=1:ncomps
            aij(i,j)=sqrt(a(i)*a(j))*(1-K(i,j));
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
        Zmidmix=Zmixture(2);
        delGRTlmix=(Zlmix-1)-log(Zlmix-Bmix)-(Amix/(2*sqrt(2)*Bmix))*log((Zlmix+(1+sqrt(2))*Bmix)/(Zlmix+(1-sqrt(2))*Bmix)) ;
        delGRTvmix=(Zvmix-1)-log(Zvmix-Bmix)-(Amix/(2*sqrt(2)*Bmix))*log((Zvmix+(1+sqrt(2))*Bmix)/(Zvmix+(1-sqrt(2))*Bmix)) ;
        delGRTmidmix=(Zmidmix-1)-log(Zmidmix-Bmix)-(Amix/(2*sqrt(2)*Bmix))*log((Zmidmix+(1+sqrt(2))*Bmix)/(Zmidmix+(1-sqrt(2))*Bmix)) ;
        if delGRTlmix<delGRTvmix && delGRTlmix <delGRTmidmix
            Zmix=Zlmix;
        elseif delGRTvmix<delGRTlmix && delGRTvmix <delGRTmidmix
            Zmix=Zvmix;
        elseif delGRTmidmix<delGRTlmix && delGRTmidmix <delGRTvmix
            Zmix=Zmidmix;
        end
    else
        Zmix=Zmixture;
    end
    
    
    %Fugacity Equation for each species in mixture
    sumterm=0;
    for i=1:ncomps
        sumterm=0;
        for j=1:ncomps
            sumterm= sumterm + y(j)*Aij(i,j);
        end
        lnphi(i)=((B(i)/Bmix)*(Zmix-1))-log(Zmix-Bmix)-(((2*sumterm)/Amix)-(B(i)/Bmix))*(Amix/(2*sqrt(2)*Bmix))*log((Zmix+(1+sqrt(2))*Bmix)/(Zmix+(1-sqrt(2))*Bmix));
        phi(i)=exp(lnphi(i));
        %         lnphi(i)=((B(i)/Bmix)*(Zvmix-1))-log(Zvmix-Bmix)-(((2*sumterm)/Amix)-(B(i)/Bmix))*(Amix/(2*sqrt(2)*Bmix))*log((Zvmix+(1+sqrt(2))*Bmix)/(Zvmix+(1-sqrt(2))*Bmix));
        %         phi(i)=exp(lnphi(i));
        f(i)=phi(i)*y(i)*P/6894.76;
    end
    
    for i=1:ncomps
        %Cordano Cubic Roots for Peng-Robinson
        [~,coeff]=cordano(A(i),B(i),Z);
        %Cubic Roots for Compressibility factor
        [Z,nbr_roots]=cubic_roots(coeff);
        if nbr_roots>1
            Zl=min(Z);
            Zv=max(Z);
            Zmid=Z(2);
            delGRTl=(Zl-1)-log(Zl-B(i))-(A(i)/(2*sqrt(2)*B(i)))*log((Zl+(1+sqrt(2))*B(i))/(Zl+(1-sqrt(2))*B(i))) ;
            delGRTv=(Zv-1)-log(Zv-B(i))-(A(i)/(2*sqrt(2)*B(i)))*log((Zv+(1+sqrt(2))*B(i))/(Zv+(1-sqrt(2))*B(i))) ;
            delGRTmid=(Zmid-1)-log(Zvmid-B(i))-(A(i)/(2*sqrt(2)*B(i)))*log((Zmid+(1+sqrt(2))*B(i))/(Zmid+(1-sqrt(2))*B(i))) ;
            
            if delGRTl<delGRTv && delGRTl <delGRTmid
                Zpure=Zl;
            elseif delGRTv<delGRTl && delGRTv <delGRTmid
                Z=Zv;
            elseif delGRTmid<delGRTl && delGRTmid <delGRTv
                Z=Zmid;
            end
        else
            Zpure=Z;
        end
        lnphi=(Z-1)-log(Z-B(i))-(A(i)/(2*sqrt(2)*B(i)))*log((Z+(1+sqrt(2))*B(i))/(Z+(1-sqrt(2))*B(i)));
        phiPure(i)=exp(lnphi);
    end
    
    %Gibbs Free Energy Minimisation
    delGmixRT=0;
    for i =1:ncomps
        phiPureMAT(counter,i)=phiPure(i);
        phiMAT(counter,i)=phi(i);
        delGmixRT=delGmixRT + y(i)*log(y(i)*phi(i)/phiPure(i));
    end
    delGmixRTMAT(counter)= delGmixRT;

    
    %DR function
    DR=0;
    for i =1:ncomps
        DR=DR + y(i)*(log(y(i)*phi(i))-log(Zfixed(i)*phiZfixed(i)));
    end
    DRMAT(counter)=DR;
%     if DR <= 1e-6
%         display(['stationary point = ' num2str(y(1))]);
%     end
    
end


%Plottin
figure(1)
plot(xC1,delGmixRTMAT,'.')
%ylim([-50,30])
xlim([0,1])
% legend('physical critical properties','optimized critical properties','physical critical properties with volume shift',0)
grid on
grid minor
% title('PV Isotherm, T=383.15K deg C')
xlabel('x1(fraction)')
ylabel('delGmixRT')

% figure(4)
% plot(xC1,DRMAT,'.','LineWidth',2)
% %ylim([-50,30])
% xlim([0,1])
% % legend('physical critical properties','optimized critical properties','physical critical properties with volume shift',0)
% grid on
% grid minor
% % title('PV Isotherm, T=383.15K deg C')
% xlabel('x1(fraction)')
% ylabel('DR function')



