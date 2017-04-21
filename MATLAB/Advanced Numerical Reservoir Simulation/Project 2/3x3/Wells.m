function [well,Qw,Qo,Q,J,Jwater] = Wells(reservoir,well,Qw,Qo,Sw,SC)

Location = well.Location;
Nx = reservoir.Nx;
Ny = reservoir.Ny;
dx = reservoir.dx;
dy = reservoir.dy;

iwell = ceil(Location(:,1)/dx);
jwell = ceil(Location(:,2)/dy);
lw = iwell + (jwell - 1) * Nx;             % l indices from x, y coordinate

Orientation = well.Orientation;            % 1=Horizontal Well,  0=Vertical Well
S = well.S;                                % skin factor
rw = well.rw;                              % well radius, ft

% Current well schedule number
if SC==1
    well.nS = well.nS + 1;
end

Type = well.Type(:,well.nS);                         % Type 1= constant rate, Type 0= constant BHP
Rate = well.Rate(:,well.nS);                         % Well Rate, ft3/d
% wRate = well.wRate(:,well.nS);                     % Well Rate, ft3/d
% oRate = well.oRate(:,well.nS);                     % Well Rate, ft3/d
BHP = well.BHP(:,well.nS);                           % BHP, psi
Wflw = well.Wflw;                                    % Previous schedule well rates (water)
Wflo = well.Wflo;                                    % Oil


% Remove effect of previous well schedule on Q (if initially Wfl is zero, no effect)
for i=1:length(lw)
    Qw(lw(i)) = Qw(lw(i)) - Wflw(i);
    Qo(lw(i)) = Qo(lw(i)) - Wflo(i);
end

kx = reservoir.kx;
ky = reservoir.ky;
kz = reservoir.kz;
h = reservoir.h;
muw = reservoir.muw;                                % constant (scalar)
muo = reservoir.muo;                                % vector
deltax = reservoir.deltax;
deltay = reservoir.deltay;

% Preallocate
Joil = sparse(Nx*Ny,Nx*Ny);
Jwater = sparse(Nx*Ny,Nx*Ny);
req = zeros(size(S));
Jw = zeros(size(S));
Jo = zeros(size(S));
Wflw = zeros(length(Type),1);
Wflo = zeros(length(Type),1);
fw = zeros(length(Type),1);
fo = zeros(length(Type),1);


for i = 1:length(Type)
    
    [krw,kro] = relperm(Sw(lw(i)));
    fw(i) = 1/(1+(kro*muw)/(krw*muo(lw(i))));
    fo(i) = 1-fw(i);
    
    if Orientation(i)==0                   % Well is vertical
        req(i) = 0.28*sqrt((ky(lw(i))/kx(lw(i)))^0.5*(deltax(lw(i)))^2 ...
            +(kx(lw(i))/ky(lw(i)))^0.5*(deltay(lw(i)))^2) / ...
            ((ky(lw(i))/kx(lw(i)))^0.25+(kx(lw(i))/ky(lw(i)))^0.25);
        

        % equivalent radius, ft
        Jw(i) = 2*pi*krw*sqrt(kx(lw(i))*ky(lw(i)))*h(lw(i))/(muw*(log(req(i)/rw(i))+S(i)))*6.33e-3;
        Jo(i) = 2*pi*kro*sqrt(kx(lw(i))*ky(lw(i)))*h(lw(i))/(muo(lw(i))*(log(req(i)/rw(i))+S(i)))*6.33e-3;
        
    else                                   % Well is horizontal
        req(i) = 0.28*sqrt((ky(lw(i))/kz(lw(i)))^0.5*(h(lw(i)))^2 ...
            +(kz(lw(i))/ky(lw(i)))^0.5*(deltay(lw(i)))^2) / ...
            ((ky(lw(i))/kz(lw(i)))^0.25+(kz(lw(i))/ky(lw(i)))^0.25);
        % equivalent radius, ft
        Jw(i) = 2*pi*krw*sqrt(ky(lw(i))*kz(lw(i)))*deltax(lw(i))/(muw*(log(req(i)/rw(i))+S(i)))*6.33e-3;
        Jo(i) = 2*pi*kro*sqrt(ky(lw(i))*kz(lw(i)))*deltax(lw(i))/(muo(lw(i))*(log(req(i)/rw(i))+S(i)))*6.33e-3;
    end

    if Type(i)==0  % well is a constant BHP well
        Jwater(lw(i),lw(i)) = Jw(i);
        Joil(lw(i),lw(i)) = Jo(i);
        Wflw(i) = Jw(i)*BHP(i);
        Wflo(i) = Jo(i)*BHP(i);
    else           % well is a constant rate well
        if Rate(i)>0   % check if its an injector
            Wflw(i) = Rate(i);
            Wflo(i) = 0;
        elseif Rate(i)<0
            Wflw(i) = fw(i)*Rate(i);
            Wflo(i) = fo(i)*Rate(i);
        end
%         Wflw(i) = wRate(i);
%         Wflo(i) = oRate(i);
    end
    
    Qw(lw(i)) = Qw(lw(i)) + Wflw(i);
    Qo(lw(i)) = Qo(lw(i)) + Wflo(i);
end

d22 = reservoir.d22;
d12 = reservoir.d12;
Q = (-d22*d12^(-1))*Qw + Qo;
J = (-d22*d12^(-1))*Jwater + Joil;

% Output in structure form
well.lw = lw;
well.Jw = Jw;
well.req = req;
well.Wflw = Wflw;
well.Wflo = Wflo;
well.fw = fw;
well.fo = fo;
well.J = J;
well.Jwater=Jwater;
end

