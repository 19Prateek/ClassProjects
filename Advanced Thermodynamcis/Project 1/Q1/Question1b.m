clear;
clc;
%Properties
P= 1e5;          %Pa
Temp = 313.15;   %K
Tc= 370;         %K
Pc= 4.25e+6;     %Pa
w= 0.152;
Rgas= 8.3144598;  %Pa m^3 K^-1 mol^-1
tol =1e-6;
Z=[0,0,0];
lnphiV=1;
lnphiL=0;
%**************************************************************************
while abs(lnphiL-lnphiV) > tol
    %Peng Robinson Equation of State
    [a,b,A,B,k,alpha]=PR_para(Pc,Tc,w,P,Temp,Rgas);
    %Cordano Cubic Roots for Peng-Robinson
    [Z,coeff]=cordano(A,B,Z);
    %Cubic Roots for Compressibility factor
    Zl=min(Z);
    Zv=max(Z);
    %Fugacity Equation
    lnphiV=(Zv-1)-log(Zv-B)-(A/(2*sqrt(2)*B))*log((Zv+(1+sqrt(2))*B)/(Zv+(1-sqrt(2))*B));
    lnphiL=(Zl-1)-log(Zl-B)-(A/(2*sqrt(2)*B))*log((Zl+(1+sqrt(2))*B)/(Zl+(1-sqrt(2))*B));
    phiV=exp(lnphiV);
    phiL=exp(lnphiL);
    P=abs(P*(phiL/phiV));
end
Vl=Zl*Rgas*Temp/P;
Vv=Zv*Rgas*Temp/P;

%P-R equation, solving for P at different Vs
V=0.000091:0.000001:Vv+4.5e-3;
aMat=a*ones(size(V));
bMat=b*ones(size(V));
TempMat=Temp*ones(size(V));
PisoT=((Rgas*TempMat)./(V-bMat))-(aMat./(V.*(V+bMat)+bMat.*(V-bMat)));
PisoTcorr=PisoT;
for i=1:size(V,2)
    if V(i)<=Vv && V(i)>= Vl
        PisoTcorr(i)=P;
    end
end

%Plotting 
disp (['Vapour Pressure = ' num2str(P/1e6) ' Mpa']);
disp([' Molar volume of liquid equilibrium phase = ' num2str(Vl) ' m^3/mol'])
disp([' Molar volume of vapour equilibrium phase = ' num2str(Vv) ' m^3/mol'])
figure(2)
semilogx(V,PisoTcorr/1e6,'LineWidth',2)
hold on
semilogx(V,PisoT/1e6,'-.','LineWidth',0.1)
ylim([1,2.2])
xlim([6e-5,3e-3])
legend('PV Isotherm T=40 deg C','Uncorrected',0)
grid on
grid minor
title('PV Isotherm, T=40 deg C')
xlabel('Molar Volume(m^3/mol)')
ylabel('Pressure (MPa)')

