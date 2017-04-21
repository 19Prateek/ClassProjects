
function [fprime,fw,Sw,xd1,xd2,xd3] = fractionalflow
muo=1.03;                                     % oil viscosity (cp)
muw=0.383211;                                 % water viscosity(cp)
Sw=(linspace(0.2,0.6,10000))';
Sw0 = 0.2;
[krw,kro]=relpermBL(Sw);                      
m=(krw*muo)./(muw*kro);                        % mobility ratio calculation
fw= 1./(1+(1./m));                             % fractional flow for water
[krw0,kro0]=relpermBL(Sw0);
m0=(krw0*muo)./(muw*kro0);                    % end point mobility ratio calculation
fw0= 1./(1+(1./m0));                          % end point fractional flow for water
% Solve for Swf from guess value read off the chart
Swf0 = 0.479;
Swf = fsolve(@myfun,Swf0);
    function F = myfun(Swf)
        dSw = 0.00001;
        Sw1 = Swf;
        Sw2 = Swf+dSw;
        [krw1,kro1] = relpermBL(Sw1);
        [krw2,kro2] = relpermBL(Sw2);
        M1 = (krw1/muw)/(kro1/muo);
        M2 = (krw2/muw)/(kro2/muo);
        fw1 = 1/(1+1/M1);
        fw2 = 1/(1+1/M2); 
        F = (fw2-fw1)/dSw - (fw1-fw0)/(Sw1-Sw0);
    end
fprime = zeros(length(fw),1);
for i=2:length(fprime)
    fprime(i) = (fw(i) - fw(i-1))/(Sw(i) - Sw(i-1));
end
[~,ind_Swf] = min(abs(Sw-Swf));
fprime(1:ind_Swf-1) = fprime(ind_Swf);
td1 = 0.1;
xd1 = td1*fprime;
td2 = 0.2;
xd2 = td2*fprime;
tdbt = 1/fprime(ind_Swf);
xd3 = tdbt*fprime;


figure(1)
plot([1;xd1],[0.2;Sw],'-','linewidth',2.5)
hold on
plot([1;xd2],[0.2;Sw],'-','linewidth',2.5)
plot([1;xd3],[0.2;Sw],'-','linewidth',2.5)
hold off
grid on
title('SATURATION PROFILES AT DIFFERENT DIMENSIONLESS TIMES')
xlabel('x_{D}')
ylabel('S_{w}')
axis([0 1 0 1])
legend('t_{D} = 0.1 (analytical)','t_{D} = 0.2(analytical)','t_{D,bt} = 0.33 (analytical)',0)

end



