function[krw,kro]=relperm(Sw)

Sor = 0.2;
Swr=0.2;
krw0=0.2;
kro0=1;
Nw=3;
No=3;
S =( Sw -Swr)/(1-Swr-Sor);

krw = krw0*(S).^Nw;
kro = kro0*(1-S).^No;


end
