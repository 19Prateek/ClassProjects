%Substance Properties
clear;
clc;
P= 2.5e6;         %Pa
Temp = 272.039;   %K
Tc   = 304.2;       %K
Pc   = 7.376e+6;   %Pa
w    = 0.225;
Rgas= 8.3144598;  %Pa m^3 K^-1 mol^-1
tol =1e-6;
Z=[0,0,0];
lnphiV=1;
lnphiL=0;
Pref=1e5;

%Peng Robinson Equation of State
[a,b,A,B,k,alpha]=PR_para(Pc,Tc,w,P,Temp,Rgas);
%Cordano Cubic Roots for Peng-Robinson
[Z,coeff]=cordano(A,B,Z);
%Cubic Roots for Compressibility factor
Zl=min(Z);
Zv=max(Z);
Vl=Zl*Rgas*Temp/P;
Vv=Zv*Rgas*Temp/P;
%Normalised Gibbs Free Energy calculation
delGRTl=(Zl-1)-log(Zl-B)-(A/(2*sqrt(2)*B))*log((Zl+(1+sqrt(2))*B)/(Zl+(1-sqrt(2))*B)) + log(P/Pref);
delGRTv=(Zv-1)-log(Zv-B)-(A/(2*sqrt(2)*B))*log((Zv+(1+sqrt(2))*B)/(Zv+(1-sqrt(2))*B)) + log(P/Pref);


disp([' Z liquid = ' num2str(Zl)])
disp([' Z vapour = ' num2str(Zv)])
disp([num2str(delGRTl) ' delGRT liquid'])
disp([num2str(delGRTv) ' delGRT vapour'])
disp('By the basic principle of minimisation of Gibbs Free Energy, ')
if (delGRTl>delGRTv)
disp('Carbon dioxide is vapour at the given conditions, because delGRTv < delGRTl. Zv results in lower Gibbs free energy ')
else
    disp('Carbon dioxide is liquid at the given conditions delGRTl < delGRTv. Zl results in lower Gibbs free energy')
end
