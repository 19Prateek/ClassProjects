function Sw = pcap(Pc)
%To initialise the reservoir saturations bsaed on capillary pressure
Swcon = 0.2;                                    % connate Water Saturation
Pe = 3.5;                                       % Entry Pressure(psi)
lambda = 2;                                     % Pore size distibution index
Sw = Swcon + ((Pc/Pe).^(-lambda))*(1-Swcon);
end