%Peng Robinson Equation of State
clear;
clc;
P= 15;     %Pa
Temp = 350;    %K
Tc= 370;       %K
Pc= 4.25e+6;   %Pa
Rgas= 8.3144598;     %Pa m^3 K^-1 mol^-1

%Peng Robinson EOS Parameters
w=0.152;
k=0.37464+1.5422*w-0.26992*w^2;
alpha=(1+k*(1-sqrt(Temp/Tc)))^2;
b=0.07780*Rgas*Tc/Pc;                   %m^3/mol
a=0.45724*(Rgas*Tc)^2*alpha/Pc;    %Pa-m^6/mol^2
fratio=2;
Z=[0,0,0];
i=0;
while abs(fratio-1) > 0.0001
%Cordano Cubic Roots for Peng-Robinson   
B=(P*b)/(Rgas*Temp);
A=(a*P)/(Rgas^2*Temp^2);
a2= -(1-B);
a1=  (A-3*B^2-2*B);
a0= -(A*B-B^2-B^3);
Q=(3*a1-a2^2)/9;
R =(9*a2*a1-27*a0-2*a2^3)/54;
D=Q^3+R^2;
coeff=[1,a2,a1,a0];
if(D<0)
    theta=acos(R/sqrt(-Q^3));
    Z(1)=2*sqrt(-Q)*cos(theta/3)-a2/3;
    Z(2)=2*sqrt(-Q)*cos((theta+2*pi())/3)-a2/3;
    Z(3)=2*sqrt(-Q)*cos((theta+4*pi())/3)-a2/3;
else
    S=(R+sqrt(D))^(1/3);
    T=(R-sqrt(D))^(1/3);
    Z(1)=-a2/3 +(S+T);
    Z(2)=-0.5*(S+T)-(1/3)*a2 + 0.5i*sqrt(3)*(S-T);
    Z(3)=-0.5*(S+T)-(1/3)*a2 - 0.5i*sqrt(3)*(S-T);
end


total_roots=roots(coeff);
if (isreal(total_roots) == 0)
    for i=1:3
        if (isreal(total_roots(i)) ==1)
            Z = total_roots(i);
        end
    end
else 
    Z_temp = sort(Z);
    if (Z_temp(1) < 0)
        Z = Z_temp(3);% if the lowest root is less than zero
        nbr_roots = 1;% there is only one positive root
    else
        Z=Z_temp;
    end
end
        
       
%Cubic Roots for Compressibility factor
Zl=min(Z);
Zv=max(Z);

fv=P*exp((Zv-1)-log(Zv-B)-(A/(2*sqrt(2)*B))*log((Zv+(1+sqrt(2)*B))/(Zv+(1-sqrt(2)*B))));
fl=P*exp((Zl-1)-log(Zl-B)-(A/(2*sqrt(2)*B))*log((Zl+(1+sqrt(2)*B))/(Zl+(1-sqrt(2)*B))));
fratio=(fl/fv);
P=P*(fl/fv);
i=i+1;
disp(i)
end

disp (P);