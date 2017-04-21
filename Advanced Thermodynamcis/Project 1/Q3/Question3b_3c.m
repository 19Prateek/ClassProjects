close all;
clear;
%Properties
clc
Temp = 272.039;                    %K
Tc   = 304.2;                      %K
Pc   = 7.376e+6;                   %Pa
w    = 0.225;
Rgas= 8.3144598;                   %Pa m^3 K^-1 mol^-1
P= 1e1;                            %Pa
tol =1e-6;
Z=[0,0,0];
lnphiV=1;
lnphiL=0;
Pref=1e5;

% ************************************Part b*******************************
V1=0.00003:1e-8:0.01;
k=0.37464+1.5422*w-0.26992*w^2;
alpha=(1+k*(1-sqrt(Temp/Tc)))^2;
b=0.07780*Rgas*Tc/Pc;              %m^3/mol
a=0.45724*(Rgas*Tc)^2*alpha/Pc;    %Pa-m^6/mol^2
%P-R EOS
PisoT=Rgas*Temp*(1./((V1)-b))-a*(1./(V1.*(V1+b)+b*(V1-b)));
BMat=(PisoT*b)/(Rgas*Temp);
AMat=(PisoT*a)/(Rgas^2*Temp^2);
Z=(PisoT.*V1)/(Rgas*Temp);
delGRT=(Z-1)-log(Z-BMat)-(AMat./(2*sqrt(2)*BMat)).*log((Z+(1+sqrt(2)).*BMat)./(Z+(1-sqrt(2)).*BMat)) + log(PisoT/Pref);

%***********************************Part c********************************

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
V2=4e-5:0.000001:Vv+4.5e-3;
aMat=a*ones(size(V2));
bMat=b*ones(size(V2));
TempMat=Temp*ones(size(V2));
PisoT=((Rgas*TempMat)./(V2-bMat))-(aMat./(V2.*(V2+bMat)+bMat.*(V2-bMat)));
PisoTcorr=PisoT;
for i=1:size(V2,2)
    if V2(i)<=Vv && V2(i)>= Vl
        PisoTcorr(i)=P;
    end
end
delGRTl=(Zl-1)-log(Zl-B)-(A/(2*sqrt(2)*B))*log((Zl+(1+sqrt(2))*B)/(Zl+(1-sqrt(2))*B)) + log(P/Pref);
delGRTv=(Zv-1)-log(Zv-B)-(A/(2*sqrt(2)*B))*log((Zv+(1+sqrt(2))*B)/(Zv+(1-sqrt(2))*B)) + log(P/Pref);

%Plotting, Results for part b and c
hold off
disp('Part b: ')
disp (['Pvap = ' num2str(P) ' Mpa']);
disp (['Vliquid = ' num2str(Vl) ' m3/mol']);
disp (['Vvapour = ' num2str(Vv) ' m3/mol']);

%Stability Regions:
unstable=delGRT;
metastable=delGRT;
figure(1)
for i=1:size(V1,2)
    if V1(i)<=Vv && V1(i)>= Vl
        if (unstable(i+1)-unstable(i))/(V1(i+1)-V1(i))<=0
            metastable(i)=delGRT(i);
            unstable(i)=NaN;
            delGRT(i)=NaN;
        else
            unstable(i)=delGRT(i);
            metastable(i)=NaN;
            delGRT(i)=NaN;
        end
    else
        unstable(i)=NaN;
        metastable(i)=NaN;
    end
end

xmeta=min(metastable);
ymeta=max(metastable);
xun=min(unstable);
yun=max(unstable);
ind1=find((metastable==xmeta));
ind2=find((metastable==ymeta));
ind3=find((unstable==xun));
ind4=find((unstable==yun));
semilogx(V1,delGRT,'-k','LineWidth',2)

hold on
semilogx(V1,unstable,'-.','LineWidth',2)
semilogx(V1,metastable,'-r','LineWidth',2)
ylim([2,7])
legend('stable','unstable','metastable',0)
grid on
grid minor
title('delG/RT vs. molar V : Stable, unstable and metastable regions')
xlabel('Molar Volume(m^3/mol)')
ylabel('delG/(RT)')
disp(['metastable region molar volume(m^3/mol) : ' num2str(Vl) ' to ' num2str(V1(ind1)) ' and ' num2str(V1(ind2)) ' to ' num2str(Vv)])
disp(['unstable region molar volume(m^3/mol): ' num2str(V1(ind3)) ' to ' num2str(V1(ind4)) ])
disp('Part c :')
disp([num2str(delGRTl) ' delGRT liquid'])
disp([num2str(delGRTv) ' delGRT vapour'])
disp('The values of molar Gibbs free energy are equal for the two phases at equilibrium')
% 
% %PV Isotherm
% figure(2)
% semilogx(V2,PisoTcorr/1e6,'LineWidth',2)
% hold on
% semilogx(V2,PisoT/1e6,'-.','LineWidth',0.1)
% legend('PV Isotherm T=30 deg F','Uncorrected',0)
% grid on
% grid minor
% title('PV Isotherm : T=30 deg F')
% xlabel('Molar Volume(m^3/mol)')
% ylabel('Pressure (MPa)')

%Liquid Vapour Region division
liquid=delGRT;
vapour=delGRT;
middle=delGRT;
figure(3)
for i=1:size(V1,2)
    if V1(i)<=Vv && V1(i)>= Vl
        middle(i)=delGRTv;
    else
        middle(i)=NaN;
    end
    if V1(i)<=Vv 
        liquid(i)=delGRT(i);
        vapour(i)=NaN;
    end
    if V1(i)>= Vl
        vapour(i)=delGRT(i);
        liquid(i)=NaN;
    end
end
semilogx(V1,liquid,'-b','LineWidth',2)
hold on
semilogx(V1,vapour,'-r','LineWidth',2)
% semilogx(V1,middle,'-b','LineWidth',2)
ylim([2,7])
grid on
grid minor
legend('Liquid Region','Vapour Region',0)
title('delG/RT vs. molar V : Liquid and Vapour Region Division')
xlabel('Molar Volume(m^3/mol)')
ylabel('delG/(RT)')

figure(4)
for i=1:size(V1,2)
    if V1(i)<=Vv && V1(i)>= Vl
        middle(i)=delGRTv;
    else
        middle(i)=NaN;
    end
    if V1(i)<=Vv 
        liquid(i)=delGRT(i);
        vapour(i)=NaN;
    end
    if V1(i)>= Vl
        vapour(i)=delGRT(i);
        liquid(i)=NaN;
    end
end

semilogx(V1,liquid,'-b','LineWidth',2)
hold on
semilogx(V1,vapour,'-b','LineWidth',2)
semilogx(V1,middle,'-b','LineWidth',2)
ylim([2,7])
grid on
grid minor
title('delG/RT vs. molar V ')
legend('delG/RT Isotherm T=30 deg F',0)
xlabel('Molar Volume(m^3/mol)')
ylabel('delG/(RT)')
