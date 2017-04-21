function [a,b,A,B,k,alpha]= PR_para(Pc,Tc,w,P,Temp,Rgas)

%Peng Robinson EOS Parameters
if w <= 0.49
    k=0.37464+ 1.5422*w-0.26992*w^2;
else
    k=0.37964+ w*(1.48503 + w*(-0.164423+0.0166666*w));
end
alpha=(1+k*(1-sqrt(Temp/Tc)))^2;
b=0.07780*Rgas*Tc/Pc;                   %m^3/mol
a=0.45724*(Rgas*Tc)^2*alpha/Pc;         %Pa-m^6/mol^2
B=(P*b)/(Rgas*Temp);
A=(a*P)/(Rgas^2*Temp^2);
end