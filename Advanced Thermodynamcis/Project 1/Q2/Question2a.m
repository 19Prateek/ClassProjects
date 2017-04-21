%****************************physical critical properties*****************
%Properties
clear;
clc;
P= 5e5;                %Pa
Temp = 383.15;         %K
Tc   = 871.16;         %K
Pc   = 5.53e5;         %Pa
w    = 1.4678;
Rgas= 8.3144598;        %Pa m^3 K^-1 mol^-1
Z=[0,0,0];
lnphiV=1;
lnphiL=0;

%Peng Robinson Equation of State
[a,b,A,B,k,alpha]=PR_para(Pc,Tc,w,P,Temp,Rgas);
%Cordano Cubic Roots for Peng-Robinson
[~,coeff]=cordano(A,B,Z);
Z=cubic_roots(coeff);
%Cubic Roots for Compressibility factor
Zl=min(Z);
Zv=max(Z);
Vl=Zl*Rgas*Temp/P;
Vv=Zv*Rgas*Temp/P;

%P-R equation, solving for P at different Vs
V1=0.000022:0.0000001:0.01;
aMat=a*ones(size(V1));
bMat=b*ones(size(V1));
TempMat=Temp*ones(size(V1));
PisoT1=((Rgas*TempMat)./(V1-bMat))-(aMat./(V1.*(V1+bMat)+bMat.*(V1-bMat)));
PisoTcorr1=PisoT1;

%********************optimized critical properties*************************

%Properties
P= 1e5;              %Pa
Temp = 383.15;       %K
Tc   = 918.25;       %K
Pc   = 9.16e5;       %Pa
w    = 1.2184;
Rgas= 8.3144598;  %Pa m^3 K^-1 mol^-1
tol =1e-2;
Z=[0,0,0];
lnphiV=1;
lnphiL=0;

%Peng Robinson Equation of State
[a,b,A,B,k,alpha]=PR_para(Pc,Tc,w,P,Temp,Rgas);
%Cordano Cubic Roots for Peng-Robinson
[~,coeff]=cordano(A,B,Z);
Z=cubic_roots(coeff);
%Cubic Roots for Compressibility factor
Zl=min(Z);
Zv=max(Z);
Vl=Zl*Rgas*Temp/P;
Vv=Zv*Rgas*Temp/P;

%P-R equation, solving for P at different Vs
V2=0.000022:0.0000001:0.01;
aMat=a*ones(size(V2));
bMat=b*ones(size(V2));
TempMat=Temp*ones(size(V2));
PisoT2=((Rgas*TempMat)./(V2-bMat))-(aMat./(V2.*(V2+bMat)+bMat.*(V2-bMat)));


%**************physical critical properties with volume shift**************
%Properties
P= 1e5;           %bar
Temp = 383.15;    %K
Tc   = 871.16;    %K
Pc   = 5.53e5;    %Pa
w    = 1.4678;
Rgas= 8.3144598;  %Pa m^3 K^-1 mol^-1
tol =1e-3;
Z=[0,0,0];
lnphiV=1;
lnphiL=0;



%Peng Robinson Equation of State
[a,b,A,B,k,alpha]=PR_para(Pc,Tc,w,P,Temp,Rgas);
%Cordano Cubic Roots for Peng-Robinson
[~,coeff]=cordano(A,B,Z);
Z=cubic_roots(coeff);
%Cubic Roots for Compressibility factor
Zl=min(Z);
Zv=max(Z);
Vl=Zl*Rgas*Temp/P;
Vv=Zv*Rgas*Temp/P;
%Volume Shift Parameters
M= 506.972e-3;         %kg/mol;
d=2.258;
e=0.1823;
s=1-(d/((M*1e3)^e));
c=s*b;


%P-R equation, solving for P at different Vs
V3=0.000022:0.0000001:0.01;
aMat=a*ones(size(V3));
bMat=b*ones(size(V3));
C=c*ones(size(V3));
TempMat=Temp*ones(size(V3));
PisoT3=((Rgas*TempMat)./(V3+C-bMat))-(aMat./((V3+C).*((V3+C)+bMat)+bMat.*((V3+C)-bMat)));

%Plotting
close all
plot(V1,PisoT1/1e6,'.','LineWidth',2)
hold on
plot(V2,PisoT2/1e6,'.','LineWidth',2)
plot(V3,PisoT3/1e6,'.','LineWidth',2)
ylim([-50,80])
xlim([1e-4,0.01])
legend('physical critical properties','optimized critical properties','physical critical properties with volume shift',0)
grid on
grid minor
title('PV Isotherm, T=383.15K deg C')
xlabel('Molar Volume(m^3/mol)')
ylabel('Pressure (MPa)')

