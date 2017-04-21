function [T,Tw,To,D,G,reservoir] = TDGBL(numerical,reservoir,P,Sw,Pc)
%TBQG calculates and assings the T, D, Q, and G matrices based on fluid
%properties, saturations and pressures from the previous timestep.
Nx=reservoir.Nx;
Ny=reservoir.Ny;
kx=reservoir.kx;
ky=reservoir.ky;
deltax=reservoir.deltax;
deltay=reservoir.deltay;
h=reservoir.h;
z=reservoir.z;
Bw = reservoir.Bw;
BoV = reservoir.Bo;
muw = reservoir.muw;
muoV = reservoir.muo;
cr = reservoir.cr;
cw = reservoir.cw;
co = reservoir.co;
rhow = reservoir.rhow;
rhoo = reservoir.rhoo;
phi=reservoir.phi;
dt = numerical.dt;
%^^^^^^^^^^^^^^^^^^^^^^Preallocation of Matrices^^^^^^^^^^^^^^^^^^^^^^^^^^^
T = sparse(Nx*Ny,Nx*Ny);
Tw = T;
To = T;
d = zeros(Nx*Ny,4);
D = sparse(Nx*Ny,Nx*Ny);
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^T,B,G Assignment^^^^^^^^^^^^^^^^^^^^^^^^^^^
for l = 1:Nx*Ny
    if mod(l-1,Nx)~=0              % The grid block is not on left boundary
        kV = [kx(l-1);kx(l)];
        dxV = [deltax(l-1);deltax(l)];
        dyV = [deltay(l-1);deltay(l)];
        hV = [h(l-1);h(l)];
        dir = 3;
        % Upwinding Scheme
        if P(l-1)<P(l)
            Bo = BoV(l);
            muo = muoV(l);
            [krw,kro]=relpermBL(Sw(l));
        else
            Bo = BoV(l-1);
            muo = muoV(l-1);
            [krw,kro]=relpermBL(Sw(l-1));
        end
        [Tw(l,l-1),To(l,l-1)] = transBL(kV,Bw,Bo,krw,kro,muw,muo,dxV,dyV,hV,dir);
        Tw(l,l) = Tw(l,l) - Tw(l,l-1);
        To(l,l) = To(l,l) - To(l,l-1);
    end
         
    if mod(l,Nx)~=0               % The grid block is not on right boundary
        kV = [kx(l);kx(l+1)];
        dxV = [deltax(l); deltax(l+1)];
        dyV = [deltay(l); deltay(l+1)];
        hV = [h(l); h(l+1)];
        dir = 2;
        % Upwinding Scheme
        if P(l)<P(l+1)
            Bo = BoV(l+1);
            muo = muoV(l+1);
            [krw,kro]=relpermBL(Sw(l+1));
        else
            Bo = BoV(l);
            muo = muoV(l);
            [krw,kro]=relpermBL(Sw(l));
        end
        [Tw(l,l+1),To(l,l+1)] = transBL(kV,Bw,Bo,krw,kro,muw,muo,dxV,dyV,hV,dir);
        Tw(l,l) = Tw(l,l) - Tw(l,l+1);
        To(l,l) = To(l,l) - To(l,l+1);
    end
    
    if l>Nx                      % The grid block is not on bottom boundary
        kV = [ky(l-Nx);ky(l)];
        dxV = [deltax(l-Nx); deltax(l)];
        dyV = [deltay(l-Nx); deltay(l)];
        hV = [h(l-Nx); h(l)];
        dir = 4;
        % Upwinding Scheme
        if P(l-Nx)<P(l)
            Bo = BoV(l);
            muo = muoV(l);
            [krw,kro]=relpermBL(Sw(l));
        else
            Bo = BoV(l-Nx);
            muo = muoV(l-Nx);
            [krw,kro]=relpermBL(Sw(l-Nx));
        end
        [Tw(l,l-Nx),To(l,l-Nx)] = transBL(kV,Bw,Bo,krw,kro,muw,muo,dxV,dyV,hV,dir);
        Tw(l,l) = Tw(l,l) - Tw(l,l-Nx);
        To(l,l) = To(l,l) - To(l,l-Nx);
        
    end
    
    if l<=Nx*(Ny-1)                 % The grid block is not on top boundary
        kV = [ky(l);ky(l+Nx)];
        dxV = [deltax(l); deltax(l+Nx)];
        dyV = [deltay(l); deltay(l+Nx)];
        hV = [h(l); h(l+Nx)];
        dir = 1;
        % Upwinding Scheme
        if P(l)<P(l+Nx)
            Bo = BoV(l+Nx);
            muo = muoV(l+Nx);
            [krw,kro]=relpermBL(Sw(l+Nx));
        else
            Bo = BoV(l);
            muo = muoV(l);
            [krw,kro]=relpermBL(Sw(l));
        end
        [Tw(l,l+Nx),To(l,l+Nx)] = transBL(kV,Bw,Bo,krw,kro,muw,muo,dxV,dyV,hV,dir);
        Tw(l,l) = Tw(l,l) - Tw(l,l+Nx);
        To(l,l) = To(l,l) - To(l,l+Nx);
        
    end
    
    V = deltax(l)*deltay(l)*h(l);
    d(l,1) = V*Sw(l)*phi(l)*(cr+cw)/(Bw*dt);             % d11 value
    d(l,2) = V*phi(l)*(1-0)/(dt*Bw);                     % d12 (NEEDS TO BE UPDATED WHEN PC IS NOT IGNORED)
    d(l,3) = V*phi(l)*(1-Sw(l))*(cr+co)/(BoV(l)*dt);     % d21 value
    d(l,4) = V/dt*(-phi(l)/BoV(l));                      % d22 value
    
    D(l,l) = -d(l,4)/d(l,2)*d(l,1)+d(l,3);        
end
d11 = sparse(diag(d(:,1)));
d12 = sparse(diag(d(:,2)));
d21 = sparse(diag(d(:,3)));
d22 = sparse(diag(d(:,4)));
T = (-d22*d12^(-1))*Tw+To;
G = 0;
reservoir.d11 = d11;
reservoir.d12 = d12;
reservoir.d21 = d21;
reservoir.d22 = d22;


end

