function Pc = capillary(Sw)
%To determine the capillary pressures at a given saturation
Pe = 3.5;                                     % Entry Pressure(psi)
Swcon = 0.2;                                  % connate Water Saturation
lambda = 2;                                   % Pore size distibution index
Pc = Pe*((Sw - Swcon)/(1-Swcon)).^(-1/lambda);
end

