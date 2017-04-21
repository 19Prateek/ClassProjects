function [Tw,To] = transBL(kV,Bw,Bo,krw,kro,muw,muo,dxV,dyV,hV,dir)
% Indexing dir (based on NEWS): 1 - North (i,j+1), 2 - East (i+1,j), 3 - West
% (i-1,j), 4 - South (i,j-1)  
% negative dir - Current block itself (flow direction based on NEWS)
if dir == 1
% transmissibility between j+1 and j at same i
    Alp1 = dxV(2) * hV(2);
    Al = dxV(1) * hV(1);
    Tw = (krw/(muw*Bw)) * (2*kV(1)*Al*kV(2)*Alp1/(kV(1)*Al*dyV(2)+kV(2)*Alp1*dyV(1)));
    To = (kro/(muo*Bo)) * (2*kV(1)*Al*kV(2)*Alp1/(kV(1)*Al*dyV(2)+kV(2)*Alp1*dyV(1)));

elseif dir == 2
% transmissibility between i and i+1 at same j
    Alp1 = dyV(2) * hV(2);
    Al = dyV(1) * hV(1);
    Tw = (krw/(muw*Bw)) * (2*kV(1)*Al*kV(2)*Alp1/(kV(1)*Al*dxV(2)+kV(2)*Alp1*dxV(1)));
    To = (kro/(muo*Bo)) * (2*kV(1)*Al*kV(2)*Alp1/(kV(1)*Al*dxV(2)+kV(2)*Alp1*dxV(1)));

elseif dir == 3
% transmissibility between i and i-1 at same j    
    Alm1 = dyV(1) * hV(1);
    Al = dyV(2) * hV(2);
    Tw = (krw/(muw*Bw)) * (2*kV(1)*Alm1*kV(2)*Al/(kV(1)*Alm1*dxV(2)+kV(2)*Al*dxV(1)));
    To = (kro/(muo*Bo)) * (2*kV(1)*Alm1*kV(2)*Al/(kV(1)*Alm1*dxV(2)+kV(2)*Al*dxV(1)));
    
elseif dir == 4
% transmissibility between j and j-1 at same i
    Alm1 = dxV(1) * hV(1);
    Al = dxV(2) * hV(2);
    Tw = (krw/(muw*Bw)) * (2*kV(1)*Alm1*kV(2)*Al/(kV(1)*Alm1*dyV(2)+kV(2)*Al*dyV(1)));
    To = (kro/(muo*Bo)) * (2*kV(1)*Alm1*kV(2)*Al/(kV(1)*Alm1*dyV(2)+kV(2)*Al*dyV(1)));
    
end
Tw = -6.33e-3*Tw;
To = -6.33e-3*To;                    
end