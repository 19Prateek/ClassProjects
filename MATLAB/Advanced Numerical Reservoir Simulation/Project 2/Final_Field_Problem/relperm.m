function[krw,kro]=relperm(Sw)
% relperm gives relative permeabilities of oil and water for a given Sw.
Swcon = 0.2;                      % connate Water Saturation
Sorw = 0.4;                       % Residual oil saturation to water
S = (Sw-Swcon)/(1-Sorw-Swcon);    % Normalised Water Saturation
krw0=0.3;                         % End point relative permeability to water
kro0=0.8;                         % End point relative permeability to oil
Nw=2;                             % Corey exponent for water
No=2;                             % Corey exponent for oil
krw = krw0*S.^Nw;
kro = kro0*(1-S).^No;
end
