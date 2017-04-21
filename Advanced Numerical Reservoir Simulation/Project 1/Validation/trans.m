function T = trans(kV,Bw,mu,dxV,dyV,hV,dir)
%Transmissibility funtion to calculate the 4 interblock transmissibilities
%associated with each grid block

% Indexing dir (based on NEWS): 1 - North (i,j+1), 2 - East (i+1,j), 3 - West
% (i-1,j), 4 - South (i,j-1)  
% negative dir - Current block itself (flow direction based on NEWS)


if dir == 1
% transmissibility between j+1 and j at same i
    Alp1 = dxV(2) * hV(2);
    Al = dxV(1) * hV(1);
    T = (1/(mu*Bw)) * (2*kV(1)*Al*kV(2)*Alp1/(kV(1)*Al*dyV(2)+kV(2)*Alp1*dyV(1)));

elseif dir == 2
% transmissibility between i and i+1 at same j
    Alp1 = dyV(2) * hV(2);
    Al = dyV(1) * hV(1);
    T = (1/(mu*Bw)) * (2*kV(1)*Al*kV(2)*Alp1/(kV(1)*Al*dxV(2)+kV(2)*Alp1*dxV(1)));

elseif dir == 3
% transmissibility between i and i-1 at same j    
    Alm1 = dyV(1) * hV(1);
    Al = dyV(2) * hV(2);
    T = (1/(mu*Bw)) * (2*kV(1)*Alm1*kV(2)*Al/(kV(1)*Alm1*dxV(2)+kV(2)*Al*dxV(1)));
    
elseif dir == 4
% transmissibility between j and j-1 at same i
    Alm1 = dxV(1) * hV(1);
    Al = dxV(2) * hV(2);
    T = (1/(mu*Bw)) * (2*kV(1)*Alm1*kV(2)*Al/(kV(1)*Alm1*dyV(2)+kV(2)*Al*dyV(1)));
    
elseif dir == -1
% Dirichlet BC on top of block     
    Al = dxV * hV;
    T = kV*Al/(mu*Bw*dyV);
    
elseif dir == -2
% Dirichlet BC on right of block    
    Al = dyV * hV;
    T = kV*Al/(mu*Bw*dxV);

elseif dir == -3
% Dirichlet BC on left of block    
    Al = dyV * hV;
    T = kV*Al/(mu*Bw*dxV);
    
elseif dir == -4
% Dirichlet BC on bottom of block    
    Al = dxV * hV;
    T = kV*Al/(mu*Bw*dyV);
    
    
end

T = 6.33e-3*T;                    % Unit conversion
    
end