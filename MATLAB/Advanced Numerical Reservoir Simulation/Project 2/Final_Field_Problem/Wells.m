function [well,Qw,Qo,Q,J] = Wells(reservoir,well,Qw,Qo,Sw,Pc)
% Returns Q, J and well properties
%^^^^^^^^^^^^^^^^^^^^^^INITIALISE VARIABLES FOR WELL CALCULATIONS^^^^^^^^^^
Nx = reservoir.Nx;                                  
Ny = reservoir.Ny;                                  
dx = reservoir.dx;                                  
dy = reservoir.dy;                                  
kx = reservoir.kx;                                  
ky = reservoir.ky;                                 
kz = reservoir.kz;                                  
h = reservoir.h;                                  
muw = reservoir.muw;                                
muoV = reservoir.muoV;                              
Bw = reservoir.Bw;
BoV = reservoir.BoV;
deltax = reservoir.deltax;                         
deltay = reservoir.deltay;                          
Location = well.Location;                         
Orientation = well.Orientation;                
S = well.S;                                       
rw = well.rw;                                      
Type = well.Type(:,well.nS);                        
Rate = well.Rate(:,well.nS);                        
BHP = well.BHP(:,well.nS);                          
Wflw = well.Wflw;                                   
Wflo = well.Wflo;                                  

%^^^^^^^^^^^^^^^^^^^^^^^^^WELL LOCATION CALCULATIONS^^^^^^^^^^^^^^^^^^^^^^^
iwell = ceil(Location(:,1)/dx);
jwell = ceil(Location(:,2)/dy);
lw = iwell + (jwell - 1) * Nx;         

%^^^^^^^^^^^^^^^^^^^^MATRIX & VECTOR PRE-ALLOCATION^^^^^^^^^^^^^^^^^^^^^^^^
for i=1:length(lw)
    Qw(lw(i)) = Qw(lw(i)) - Wflw(i);
    Qo(lw(i)) = Qo(lw(i)) - Wflo(i);
end

Joil = sparse(Nx*Ny,Nx*Ny);
Jwater = sparse(Nx*Ny,Nx*Ny);
req = zeros(size(S));
Jw = zeros(size(S));
Jo = zeros(size(S));
Wflw = zeros(length(Type),1);
Wflo = zeros(length(Type),1);
fw = zeros(length(Type),1);
fo = zeros(length(Type),1);


%^^^^^^^^^^^^^WELL INDEX AND EQUIVALENT RADIUS CALCULATION^^^^^^^^^^^^^^^^^
for i = 1:length(Type)
    muo = muoV(lw(i));
    Bo = BoV(lw(i));
    [krw,kro] = relperm(Sw(lw(i)));
    fw(i) = 1/(1+(kro*muw)/(krw*muo));
    fo(i) = 1-fw(i);
    
    if Orientation(i)==0  % Req(ft) and J(ft^3/day-psi) for vertical wells 
        
        req(i) = 0.28*sqrt((ky(lw(i))/kx(lw(i)))^0.5*(deltax(lw(i)))^2+(kx(lw(i))/ky(lw(i)))^0.5*(deltay(lw(i)))^2)/((ky(lw(i))/kx(lw(i)))^0.25+(kx(lw(i))/ky(lw(i)))^0.25);  
        Jw(i) = 2*pi*krw*sqrt(kx(lw(i))*ky(lw(i)))*h(lw(i))/(Bw*muw*(log(req(i)/rw(i))+S(i)))*6.33e-3;
        Jo(i) = 2*pi*kro*sqrt(kx(lw(i))*ky(lw(i)))*h(lw(i))/(Bo*muo*(log(req(i)/rw(i))+S(i)))*6.33e-3;
                                                 
        
    else    % Req(ft and J(ft^3/day-psi) for Horizontal wells
        
        req(i) = 0.28*sqrt((ky(lw(i))/kz(lw(i)))^0.5*(h(lw(i)))^2 +(kz(lw(i))/ky(lw(i)))^0.5*(deltay(lw(i)))^2)/((ky(lw(i))/kz(lw(i)))^0.25+(kz(lw(i))/ky(lw(i)))^0.25);
        Jw(i) = 2*pi*krw*sqrt(ky(lw(i))*kz(lw(i)))*deltax(lw(i))/(Bw*muw*(log(req(i)/rw(i))+S(i)))*6.33e-3;
        Jo(i) = 2*pi*kro*sqrt(ky(lw(i))*kz(lw(i)))*deltax(lw(i))/(Bo*muo*(log(req(i)/rw(i))+S(i)))*6.33e-3;
                                                    
    end

    if Type(i)==0                                       % Constant BHP well                             
        Jwater(lw(i),lw(i)) = Jw(i);
        Joil(lw(i),lw(i)) = Jo(i);
        Wflw(i) = Jw(i)*BHP(i);
        Wflo(i) = Jo(i)*BHP(i);
        
    else                                               % Constant rate well                                       
        
        if Rate(i)>0                                   % Water injector                              
            Wflw(i) = Rate(i);
            Wflo(i) = 0;
        elseif Rate(i)<0                               % Producer, Rates divided by fractional flow                        
            Wflw(i) = fw(i)*Rate(i);
            Wflo(i) = fo(i)*Rate(i);
        end
        
    end
    
    Qw(lw(i)) = Qw(lw(i)) + Wflw(i);
    Qo(lw(i)) = Qo(lw(i)) + Wflo(i);
end

d22 = reservoir.d22;
d12inv = reservoir.d12inv;
well.lw = lw;
well.Jw = Jw;
well.Jo = Jo;
well.Jwater = Jwater;
well.Joil = Joil;
well.req = req;
well.Wflw = Wflw;
well.Wflo = Wflo;
well.fw = fw;
well.fo = fo;

%^^^^^^^^^^^^^^^^^^^Q(VECTOR),J(MATRIX) ASSIGNMENT^^^^^^^^^^^^^^^^^^^^^^^^^

Q = (-d22*d12inv)*Qw-(d22*d12inv)*Jwater*Pc + Qo;
J = (-d22*d12inv)*Jwater + Joil;

end

