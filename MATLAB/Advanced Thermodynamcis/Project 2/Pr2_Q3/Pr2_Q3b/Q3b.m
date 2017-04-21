clear;
% close all;
clc;
%Properties
close all


gammaL=[0 0 0];
delgrtL=gammaL;
counterL=[0 1];
v=1;
%Simulation Properties
Temp = 355.372;                  %K
P=  1.1721e+7;                   %Pa
Tc= [304.2 425.12 617.70 ];       %K
Pc= [ 73.765e5 30.70e5 21.10e5  ];     %Pa
w= [ 0.225 0.2511 0.4898 ];
Rgas= 8.3144598;  %Pa m^3 K^-1 mol^-1
ncomps=3;
%Binary Interaction Parameters
K=[0 0.12 0.1141;
   0.12 0 0
   0.1141 0 0];
aij=K;
Z=[0,0,0];        %Cubic Roots Initialisation
lnphi=Tc;
f=Tc;
phi=Tc;
%Fixed Composition Initialisation

Z=[0,0,0];
% Component 1:CO2, Component 2:nC4, Component 3:nC7, Component 4:nC
Zi=[0 0 0];
aij=K;

%**********************************
%Peng Robinson EOS parameters
[a,b,A,B,k,alpha]=PR_para(Pc,Tc,w,P,Temp,Rgas,ncomps);

comp1=[0.94067    0.059327           0];
comp2=[0     0.33139     0.66861];
vComp=[0.9106 0.0680 0.0214];
lComp=[0.7247 0.1218 0.1535];

%**************************************************************************
gamma=0:0.0001:1;
delGRTMAT=gamma;

delGRT=0;
for counter=1:(size(gamma,2))
    Zi=gamma(counter)*comp1 +(1-gamma(counter))*comp2;
    
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
        if delGRTlmix<delGRTvmix
            Zmix=Zlmix;          
        elseif delGRTvmix<delGRTlmix 
            Zmix=Zvmix;             
        end
    else
        Zmix=Zmixture;
    end  
        sumterm=0;
    for i=1:ncomps
        sumterm=0;
        for j=1:ncomps
            sumterm= sumterm + y(j)*Aij(i,j);
        end
        lnphi(i)=((B(i)/Bmix)*(Zmix-1))-log(Zmix-Bmix)-(((2*sumterm)/Amix)-(B(i)/Bmix))*(Amix/(2*sqrt(2)*Bmix))*log((Zmix+(1+sqrt(2))*Bmix)/(Zmix+(1-sqrt(2))*Bmix));
        phi(i)=exp(lnphi(i));
    end
    
        delGRT=0;
    for i =1:ncomps
        delGRT=delGRT + y(i)*log(y(i)*phi(i));
    end
    %Gibbs Free Energy Minimisation
    delGRTMAT(counter)= delGRT;
    if norm(abs(vComp-Zi)) < 0.8e-4 || norm(abs(lComp-Zi))< 1e-4
        disp([' equilibrium conc at gamma = ' num2str(gamma(counter))])
        disp([' Zi = ' num2str(Zi)]);
        delgrtL(v)=exp(delGRT);
        gammaL(v)=gamma(counter);
        counterL(v)=counter;
        v=v+1;
    end


end
g=0;
 derivative=gamma;
for g=1:(size(gamma,2)-1)
derivative(g+1)=(delGRTMAT(g+1)-delGRTMAT(g))/(gamma(g+1)-gamma(g));
end
disp(['first derivative at gamma 0.770 = ' num2str(derivative(counterL(1)))]);
disp(['first derivative at gamma 0.968 = ' num2str(derivative(counterL(2)))]);

figure(1)
%Plottin
plot(gamma,delGRTMAT,'.')
hold on
%ylim([-50,30])
xlim([0,1])

grid on
grid minor
xlabel('gamma(mixing ratio)')
ylabel('G/RT')


%Plottin
figure (2)
plot(gamma,derivative,'.')
xlim([0,1])
grid on
grid minor
xlabel('gamma(mixing ratio)')
ylabel('first derivative of G/RT')



