function [T,Tw,To,D,G,reservoir] = TDG(numerical,reservoir,P,Sw,Pc)
%TBG calculates and assings the T, D and G matrices based on fluid
%properties, saturations and pressures from the previous timestep.

%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^READ FROM USER INPUTS^^^^^^^^^^^^^^^^^^^^^
dt = numerical.dt;
Nx = reservoir.Nx;
Ny = reservoir.Ny;
kx = reservoir.kx;
ky = reservoir.ky;
phi = reservoir.phi;
deltax = reservoir.deltax;
deltay = reservoir.deltay;
h = reservoir.h;
z = reservoir.z;
Bw = reservoir.Bw;
BoV = reservoir.BoV;
muw = reservoir.muw;
muoV = reservoir.muoV;
cr = reservoir.cr;
cw = reservoir.cw;
coV = reservoir.coV;
rhow = reservoir.rhow;
rhoo = reservoir.rhoo;
T = sparse(Nx*Ny,Nx*Ny);
Tw = T;
To = T;
d = zeros(Nx*Ny,5);
D = sparse(Nx*Ny,Nx*Ny);

%^^^^^^^^^^^^^^^^^^^^^^^^^^^^T,D,G Matrix Assignment^^^^^^^^^^^^^^^^^^^^^^^
for l = 1:Nx*Ny
    % The grid block is not on left boundary
    if mod(l-1,Nx)~=0                               
        kV = [kx(l-1);kx(l)];
        dxV = [deltax(l-1);deltax(l)];
        dyV = [deltay(l-1);deltay(l)];
        hV = [h(l-1);h(l)];
        dir = 3;
        % Upwinding of relative mobility
        if P(l-1)<P(l)
            Bo = BoV(l);
            muo = muoV(l);
            [krw,kro]=relperm(Sw(l));
        else
            Bo = BoV(l-1);
            muo = muoV(l-1);
            [krw,kro]=relperm(Sw(l-1));
         end
        [Tw(l,l-1),To(l,l-1)] = trans(kV,Bw,Bo,krw,kro,muw,muo,dxV,dyV,hV,dir);
        Tw(l,l) = Tw(l,l) - Tw(l,l-1);
        To(l,l) = To(l,l) - To(l,l-1);
        
    end
    
    % The grid block is not on right boundary
    if mod(l,Nx)~=0                                
        kV = [kx(l);kx(l+1)];
        dxV = [deltax(l); deltax(l+1)];
        dyV = [deltay(l); deltay(l+1)];
        hV = [h(l); h(l+1)];
        dir = 2;
        % Upwinding of relative mobility
        if P(l)<P(l+1)
            Bo = BoV(l+1);
            muo = muoV(l+1);
            [krw,kro]=relperm(Sw(l+1));
        else
            Bo = BoV(l);
            muo = muoV(l);
            [krw,kro]=relperm(Sw(l));
        end
        [Tw(l,l+1),To(l,l+1)] = trans(kV,Bw,Bo,krw,kro,muw,muo,dxV,dyV,hV,dir);
        Tw(l,l) = Tw(l,l) - Tw(l,l+1);
        To(l,l) = To(l,l) - To(l,l+1);

    end
    
    % The grid block is not on bottom boundary
    if l>Nx                                         
        kV = [ky(l-Nx);ky(l)];
        dxV = [deltax(l-Nx); deltax(l)];
        dyV = [deltay(l-Nx); deltay(l)];
        hV = [h(l-Nx); h(l)];
        dir = 4;
        % Upwinding of relative mobility
        if P(l-Nx)<P(l)
            Bo = BoV(l);
            muo = muoV(l);
            [krw,kro]=relperm(Sw(l));
        else
            Bo = BoV(l-Nx);
            muo = muoV(l-Nx);
            [krw,kro]=relperm(Sw(l-Nx));
        end        
        [Tw(l,l-Nx),To(l,l-Nx)] = trans(kV,Bw,Bo,krw,kro,muw,muo,dxV,dyV,hV,dir);
        Tw(l,l) = Tw(l,l) - Tw(l,l-Nx);
        To(l,l) = To(l,l) - To(l,l-Nx);
        
    end
    
    % The grid block is not on top boundary
    if l<=Nx*(Ny-1)                                 
        kV = [ky(l);ky(l+Nx)];
        dxV = [deltax(l); deltax(l+Nx)];
        dyV = [deltay(l); deltay(l+Nx)];
        hV = [h(l); h(l+Nx)];
        dir = 1;
        % Upwinding of relative mobility
        if P(l)<P(l+Nx)
            Bo = BoV(l+Nx);
            muo = muoV(l+Nx);
            [krw,kro]=relperm(Sw(l+Nx));
        else
            Bo = BoV(l);
            muo = muoV(l);
            [krw,kro]=relperm(Sw(l));
        end              
        [Tw(l,l+Nx),To(l,l+Nx)] = trans(kV,Bw,Bo,krw,kro,muw,muo,dxV,dyV,hV,dir);
        Tw(l,l) = Tw(l,l) - Tw(l,l+Nx);
        To(l,l) = To(l,l) - To(l,l+Nx);
    
    end
    
    % Numerical Derivative of Pc
    dSw = 0.00001;
    Pcprime = (capillary(Sw(l)+dSw) - capillary(Sw(l)))/dSw;
    if Pcprime >= 1000000
        Pcprime = 1000000;
    end
 
    V = deltax(l)*deltay(l)*h(l);                             
    d(l,1) = V*Sw(l)*phi(l)*(cr+cw)/(Bw*dt);                % d11 value
    d(l,2) = V*phi(l)*(1-Sw(l)*phi(l)*cw*Pcprime)/(dt*Bw);  % d12 value
    d(l,3) = V*phi(l)*(1-Sw(l))*(cr+coV(l))/(BoV(l)*dt);    % d21 value
    d(l,4) = V/dt*(-phi(l)/BoV(l));                         % d22 value
    d(l,5) = 1/d(l,2);                                      % d12^-1 value
    
    % Set inverse of d12 to 0 if it's a division by 0
    if d(l,2) == 0
        d(l,5) = 0;                                         % d12^-1 value
    else
        d(l,5) = 1/d(l,2);                                  % d12^-1 value
    end
    D(l,l) = -d(l,4)/d(l,2)*d(l,1)+d(l,3);          % D matrix value
    
end

% Convert to sparse diagonal matrices for use in equations
d11 = sparse(diag(d(:,1)));
d12 = sparse(diag(d(:,2)));
d21 = sparse(diag(d(:,3)));
d22 = sparse(diag(d(:,4)));
d12inv = sparse(diag(d(:,5)));
reservoir.d11 = d11;
reservoir.d12 = d12;
reservoir.d12inv = d12inv;
reservoir.d21 = d21;
reservoir.d22 = d22;

%^^^^^^^^^^^^^^^^^^^^^^^^^^^^ACCOUNTING FOR INACTIVE BLOCKS^^^^^^^^^^^^^^
NanD_ind = isnan(D);
D(NanD_ind) = 1;                                    % Set NaN to 1 in D
NanTw_ind = isnan(Tw);
Tw(NanTw_ind) = 0;                                  % Set NaN to 0 in Tw
NanTo_ind = isnan(To);
To(NanTo_ind) = 0;                                  % Set NaN to 0 in To

T = (-d22*d12inv)*Tw+To;
G = (-d22*d12inv)*Tw*Pc+(-d22*d12inv*rhow*Tw+rhoo*To)*z;

end

