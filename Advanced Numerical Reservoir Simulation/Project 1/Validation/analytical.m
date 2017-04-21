function [qd,td,q,t] = analytical(dt,t)

% To get analytical value of rates and plot dimensionless rates vs.
% dimensionless time

phi = 0.20;        % porosity (fraction)
k = 10;            % permeability (mD)
mu = 1;            % viscosity (cp)                                                                                                                                                                                                                                                                                                               
ct = 1e-5;         % compressibilty (1/psi)
xe  = 2550;        % drainiage boundary(ft)
xf = 50;           % Half length of the Fracture(ft)
xwf = 0;           % wellbore x coordinate(ft)
Pwf = 1000;        % BHP (psi) 
Pi =  3000;        % Initial reservoir pressure(psi)
L= 2*(xe);         % Distance between two fractures(ft)
h= 50;             % Thickness of the reservoir(ft)


t=(0:dt:t)';
t=t(2:end);
td=zeros(size(t));
qd=zeros(size(t));
q=zeros(size(t));

for i=1:length(t)
    
    td(i) = 0.006328*k*t(i)/(phi*mu*ct*(xwf-xe)^2);
    qd(i) = 2*exp(-td(i)*pi^2/4)+erfc(3/2*pi*sqrt(td(i)))/sqrt(pi*td(i));
    q(i) = (k*h*(Pi-Pwf)*xf)*qd(i)/(158.0208*mu*L);
    
end



%=================Plotting dimensionless rates and time====================
% figure(1)
% loglog(td,qd,'r*')
% title('Production Profile (Dimensionless)')
% xlabel('t_{D}')
% ylabel('q_{D}')
% grid on
% 
% figure(2)
% loglog(t,q/5.615,'r*')
% title('Production Profile')
% xlabel('t (days)')
% ylabel('q (bbl/day)')
% grid on


end

