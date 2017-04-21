% 2-D Reservoir Simulation Project 1
% PGE 392K: Numerical Simulation of Reservoirs

%========SIMULATION AND RESERVOIR PARAMETERS,HETEROGENEITY & GRAVITY======
dt = 0.01;                           % Time-step size, day
t0 = 0;                              % Initial time (when initial pressure conditions apply), days
t_fin = 2052;                        % final time to evaluate pressures, days
nt = round((t_fin-t0)/dt);           % Number of time steps
rho = 0;                             % Density of water w/ gravity, psi/ft
Nx = 51;
Ny = 1;
dx = 50;                             % Block width in x-dir, ft
dy = 50;                             % Block width in y-dir, ft
deltax = dx*ones(Nx*Ny,1);           % dx Vector for use in trans
deltay = dy*ones(Nx*Ny,1);           % dy Vector for use in trans
h = 50*ones(Nx*Ny,1);
kx = zeros(Nx*Ny,1);
kx(1) = 5000;
kx(2:end)=10;
ky=kx;
phi = 0.20*ones(Nx*Ny,1);
mu = 1;                             % viscosity, cP
ct = 1e-5*ones(Nx*Ny,1);            % total compressibility, psi^-1
Bw = 1;                             % Formation volume factor
P0=3000*ones(Nx*Ny,1);
z=zeros(Nx*Ny,1);

%========Preallocate matrices and vectors==================================
T = sparse(Nx*Ny,Nx*Ny);
B = sparse(Nx*Ny,Nx*Ny);
P = zeros(Nx*Ny,nt);                  % P matrix for recording pressures
t = zeros(nt,1);                      % t vector for recording times
%===============================WELLS=====================================

Location = [25 25];
Type = [0]';               % Type 1= constant rate, type0= constant BHP
Rate = [0]';               % Well Rate, ft3/d
BHP = [1000]';             % BHP, psi
S = [0]';                  % skin factor
rw = [0.25]';              % well radius, ft
iwell = ceil(Location(:,1)/dx);
jwell = ceil(Location(:,2)/dy);
lw = iwell + (jwell - 1) * Nx;     % l indices from x, y coordinate
req = zeros(size(S));
Jw = zeros(size(S));

%Pre-allocate Q, J and Wfl matrices
Q = sparse(Nx*Ny,1);
J = sparse(Nx*Ny,Nx*Ny);
Wfl = sparse(Nx*Ny,1);


for i = 1:length(Type)
    
    req(i) = 0.28*sqrt((ky(lw(i))/kx(lw(i)))^0.5*(deltax(lw(i)))^2 ...
        +(kx(lw(i))/ky(lw(i)))^0.5*(deltay(lw(i)))^2) / ...
        ((ky(lw(i))/kx(lw(i)))^0.25+(kx(lw(i))/ky(lw(i)))^0.25);
                                                  % equivalent radius, ft
    Jw(i) = 2*pi*sqrt(kx(lw(i))*ky(lw(i)))*h(lw(i))/(mu*(log(req(i)/rw(i))+S(i)))*6.33e-3;
                                                  % Productivity index for BHP well
    if Type(i)==0       
        J(lw(i),lw(i)) = Jw(i);
        Wfl(lw(i)) = J(lw(i),lw(i))*BHP(i);
        
    else
        Wfl(lw(i)) = Rate(i);
    end
    
end

%==========================BOUNDARY CONDITIONS============================
Px1 = 1000; % Left Dirichlet boundary pressure, psi
Px2 = 1200; % Right Dirichlet boundary pressure, psi
Py1 = 3000; % Bottom Dirichlet boundary pressure, psi
Py2 = 4000; % Top Dirichlet boundary pressure, psi
wx1 = 1; % Left boundary BC type (0 - Neumann, 1 - Dirichlet)
wx2 = 0; % Right boundary BC type
wy1 = 0; % Bottom boundary BC type
wy2 = 0; % Top boundary BC type

%---------------------Calculate T,B,Q,G----------------------------------
for l = 1:Nx*Ny

    if mod(l-1,Nx)~=0
    % If block L is not on left...
        
        kV = [kx(l-1);kx(l)];
        dxV = [deltax(l-1);deltax(l)];
        dyV = [deltay(l-1);deltay(l)];
        hV = [h(l-1);h(l)];
        T(l,l-1) = -trans(kV,Bw,mu,dxV,dyV,hV,3);
        T(l,l) = T(l,l) - T(l,l-1);
    
    elseif wx1 == 1
    % If left side of block is Dirichlet BC... 
    
        kV = kx(l);
        dxV = deltax(l);
        dyV = deltay(l);
        hV = h(l);
        T(l,l) = T(l,l) + 2*trans(kV,Bw,mu,dxV,dyV,hV,-3);
        Q(l)= Q(l) + 2*trans(kV,Bw,mu,dxV,dyV,hV,-3)*(Px1-rho*z(l));
        
    elseif wx1 == 0
    % If left side of block is Neumann BC...
    
    end
    
    if mod(l,Nx)~=0
    % If block L is not on right...
    
        kV = [kx(l);kx(l+1)];
        dxV = [deltax(l); deltax(l+1)];
        dyV = [deltay(l); deltay(l+1)];
        hV = [h(l); h(l+1)];
        T(l,l+1) = -trans(kV,Bw,mu,dxV,dyV,hV,2);
        T(l,l) = T(l,l) - T(l,l+1);
    
    elseif wx2 == 1
    % If right side of block is Dirichlet BC... 
    
        kV = kx(l);
        dxV = deltax(l);
        dyV = deltay(l);
        hV = h(l);
        T(l,l) = T(l,l) + 2*trans(kV,Bw,mu,dxV,dyV,hV,-2);
        Q(l)= Q(l) + 2*trans(kV,Bw,mu,dxV,dyV,hV,-2)*(Px2-rho*z(l));
        
    elseif wx2 == 0
    % If right side of block is Neumann BC...
    
    end
    
    if l>Nx
    % If block L is not on bottom...
    
        kV = [ky(l-Nx);ky(l)];
        dxV = [deltax(l-Nx); deltax(l)];
        dyV = [deltay(l-Nx); deltay(l)];
        hV = [h(l-Nx); h(l)];
        T(l,l-Nx) = -trans(kV,Bw,mu,dxV,dyV,hV,4);
        T(l,l) = T(l,l) - T(l,l-Nx);
    
    elseif wy1 == 1
    % If bottom side of block is Dirichlet BC...
    
        kV = ky(l);
        dxV = deltax(l);
        dyV = deltay(l);
        hV = h(l);
        T(l,l) = T(l,l) + 2*trans(kV,Bw,mu,dxV,dyV,hV,-4);
        Q(l)= Q(l) + 2*trans(kV,Bw,mu,dxV,dyV,hV,-4)*(Py1-rho*z(l));
        
    elseif wy1 == 0
    % If bottom side of block is Neumann BC...
    
    end
    
    if l<=Nx*(Ny-1)
    % If block L is not on top...
    
        kV = [ky(l);ky(l+Nx)];
        dxV = [deltax(l); deltax(l+Nx)];
        dyV = [deltay(l); deltay(l+Nx)];
        hV = [h(l); h(l+Nx)];
        T(l,l+Nx) = -trans(kV,Bw,mu,dxV,dyV,hV,1);
        T(l,l) = T(l,l) - T(l,l+Nx);
    
    elseif wy2 == 1
    % If top side of block is Dirichlet BC...
    
        kV = ky(l);
        dxV = deltax(l);
        dyV = deltay(l);
        hV = h(l);
        T(l,l) = T(l,l) + 2*trans(kV,Bw,mu,dxV,dyV,hV,-1);
        Q(l)= Q(l) + 2*trans(kV,Bw,mu,dxV,dyV,hV,-1)*(Py2-rho*z(l));
        
    elseif wy2 == 0
    % If top side of block is Neumann BC...
    
    end
    
    B(l,l) = deltax(l)*deltay(l)*h(l)*phi(l)*ct(l)/Bw;
    Q(l) = Q(l) + Wfl(l);
           
end

G = 0;                                   % Gravity vector
%-------------------------------------------------------------------------
q = zeros(nt,length(lw));
Pwf = zeros(nt,length(lw));
%B=zeros(Nx*Ny,Nx*Ny);                         % Only if Steady-State
A = T+J+B/dt;
B = sparse(B);
Q = sparse(Q);
A = sparse(A);
Pold = P0;                           

for count=1:nt
    
    t(count) = count*dt+t0;
    
    Pnew = mldivide(A,B/dt*Pold+Q+G);
    
    P(:,count) = Pnew;
    
    for i=1:length(lw)
    
       if Type(i)==0
          
           q(count,i) = Jw(i)*(BHP(i)-P(lw(i),count));
           Pwf(count,i) = BHP(i);
           
       else
           
           q(count,i) = Rate(i);
           Pwf(count,i) = Rate(i)/Jw(i) + P(lw(i),count);
           
       end
           
        
    end
    
    Pold = Pnew;    

end

%============================VALIDATION PARAMETERS=========================
xe  = 2550;        % drainiage boundary(ft)
xf = 50;           % Half length of the Fracture(ft)
xwf = 0;           % wellbore x coordinate(ft)
L= 2*(xe);         % Distance between two fractures(ft)
P1=(P(1,:))';
q1=Jw(1)*(P1-BHP(1));                                     
qd1=(158.0208*q1*mu*L)/(kx(2)*h(lw(1))*(P0(lw(1))-BHP(1))*xf);       
t1 = t;
td1=(0.006328*kx(2)*t1)/(phi(lw(1))*mu*ct(lw(1))*(xe-xwf)^2);

[qd,td,q,t] = analytical(0.01,2055);

figure(3)
loglog(td1,qd1,'linewidth',1.5)
hold on
grid on
loglog(td,qd,'r*','markersize',0.95)
title('Dimensionless production rate')
xlabel('t_{D}')
ylabel('q_{D}')
legend('Simulation','Analytical')
hold off

