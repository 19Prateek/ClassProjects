function[krw,kro]=relpermBL(Sw)
Swcon = 0.2;
Swcrit = Swcon;
Sorw = 0.4;
Soirw = Sorw;
krw0=0.3;
kro0=0.8;
Nw=2;
No=2;
So=1-Sw;
krw = krw0*((Sw-Swcrit)/(1-Swcrit-Soirw)).^Nw;
kro = kro0*((So-Sorw)/(1-Swcon-Sorw)).^No;
end
