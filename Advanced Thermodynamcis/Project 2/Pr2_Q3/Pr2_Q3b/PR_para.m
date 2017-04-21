function [ai,bi,Ai,Bi,k,alpha]= PR_para(Pc,Tc,w,P,Temp,Rgas,ncomps)

k=w;
for i=1:ncomps
    %Peng Robinson EOS Parameters
    if w(i) <= 0.49
        k(i)=0.37464+ 1.5422*w(i)-0.26992*w(i)^2;
    else
        k(i)=0.37964+ w(i)*(1.48503 + w(i)*(-0.164423+0.0166666*w(i)));
    end
end
alpha=(1+k.*(1-sqrt(Temp./Tc))).^2;
bi=0.07780*Rgas*Tc./Pc;                   %m^3/mol
ai=0.45724*(Rgas*Tc).^2.*alpha./Pc;         %Pa-m^6/mol^2
Bi=(P.*bi)/(Rgas*Temp);
Ai=(ai.*P)/(Rgas^2*Temp^2);
end