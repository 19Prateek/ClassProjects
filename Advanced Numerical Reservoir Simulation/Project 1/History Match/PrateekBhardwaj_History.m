function [tc,qc]= PrateekBhardwaj_History(kratio,comp)% 2-D Reservoir Simulation Project 1

% PGE 392K: Numerical Simulation of Reservoirs
% Prateek Bhardwaj
tic        

dt = 1;                                  % Time-step size, day
t0 = 0;                                  % Initial time , days
t_fin = 1000;                            % final time , days
nt = round((t_fin-t0)/dt);               % Number of time steps
rho = 0.433;                             % Pressure gradient(water), psi/ft
zM = dlmread('PJ1-Depth.txt');           % Depth Matrix
kM = dlmread('PJ1-Permeability.txt');    % Perm Matrix
phiM = dlmread('PJ1-Porosity.txt');      % Porosity Matrix
hM = dlmread('PJ1-Thickness.txt');       % Thickness Matrix
Nx=80;Ny=75;                      % Number of blocks in x and y dir
mu = 1.95;                               % viscosity, cP
ct = comp*1e-5*ones(Nx*Ny,1);                 % total compressibility, psi^-1
Bw = 1;                                  % Formation volume factor                              
                                          
                                          
%******************************MATRIX TO VECTORS***************************
z = zeros(6000,1);                     % Depth , ft
kx = zeros(6000,1);                    % Permeability in x-dir, mD
phi = zeros(6000,1);                   % Porostiy,fraction
h = zeros(6000,1);                     % Thickness, ft

for j=1:Ny
    z((1+(j-1)*Nx):(Nx+(j-1)*Nx))=(zM(j,:))';
    kx((1+(j-1)*Nx):(Nx+(j-1)*Nx))=(kM(j,:))';
    phi((1+(j-1)*Nx):(Nx+(j-1)*Nx))=(phiM(j,:))';
    h((1+(j-1)*Nx):(Nx+(j-1)*Nx))=(hM(j,:))';
end
z=-z;                                % z is positive downwards

%***************************INITIALISATION*********************************

Lx = 6000;                          % Length in x-dir of reservoir, ft
Ly = 7500;                          % Length in y-dir of reservoir, ft
dx = Lx/Nx;                         % Block width in x-dir, ft
dy = Ly/Ny;                         % Block width in y-dir, ft
deltax = dx*ones(Nx*Ny,1);          % dx Vector for use in trans
deltay = dy*ones(Nx*Ny,1);          % dy Vector for use in trans
xpos = zeros(Nx,1);                 %(x,y) of cell centres (Uniform grids)
ypos = zeros(Ny,1);
xpos(1) = dx/2;
ypos(1) = dy/2;
for i=2:Nx
    xpos(i)=(i-1)*dx+xpos(1);
end

for i=2:Ny
    ypos(i)=(i-1)*dy+ypos(1);
end
P0max = 3500;
zmax = max(z);
P0 = P0max-rho*(zmax-z);                 % Initial reservoir pressure, psi
%****************************PREALLOCATION*********************************
T = sparse(Nx*Ny,Nx*Ny);
B = sparse(Nx*Ny,Nx*Ny);
P = zeros(Nx*Ny,nt);                  % P matrix for recording pressures
t = zeros(nt,1);                      % t vector for recording times
Q = sparse(Nx*Ny,1);
 ky = kratio*kx;                                  % Permeability in y-dir, mD
 kz = kx;                                      % Permeability in y-dir, mD          
%*******************************WELLS**************************************
Location = [4184 1569; 
            4968.5 2510.4; 
            3294.9 2928.8; 
            2562.7 4393.2; 
            1307.5 2824.2;
            890 895; 
            965 895; 
            1040 895];

iwell = ceil(Location(:,1)/dx);
jwell = ceil(Location(:,2)/dy);
lw = iwell + (jwell - 1) * Nx;             % l indices from x, y coordinate


Type = [ 1 1 1 0 0 1 1 1]';                % Type 1= constant rate, Type 0= constant BHP
Orientation = [ 0 0 0 0 0 0 0 0]';         %  1=Horizontal Well,  0=Vertical Well
Rate = [-125 -175 -750 0 0 -333.33 -333.33 -333.33]';    % Well Rate, ft3/d
BHP = [ 0 0 0 1200 1200 0 0 0]';           % BHP, psi
S = [ 6 0 0 0 6 0 0 0]';                   % skin factor
rw = [0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25]';% well radius, ft


J = sparse(Nx*Ny,Nx*Ny);
Wfl = sparse(Nx*Ny,1);
req = zeros(size(S));
Jw = zeros(size(S));

for i = 1:length(Type)
    
   if Orientation(i)==0                   % Well is vertical 
    req(i) = 0.28*sqrt((ky(lw(i))/kx(lw(i)))^0.5*(deltax(lw(i)))^2 ...
        +(kx(lw(i))/ky(lw(i)))^0.5*(deltay(lw(i)))^2) / ...
        ((ky(lw(i))/kx(lw(i)))^0.25+(kx(lw(i))/ky(lw(i)))^0.25);
                                          % equivalent radius, ft
    Jw(i) = 2*pi*sqrt(kx(lw(i))*ky(lw(i)))*h(lw(i))/(mu*(log(req(i)/rw(i))+S(i)))*6.33e-3;
      
   else                                   % Well is horizontal
        req(i) = 0.28*sqrt((ky(lw(i))/kz(lw(i)))^0.5*(h(lw(i)))^2 ...
        +(kz(lw(i))/ky(lw(i)))^0.5*(deltay(lw(i)))^2) / ...
        ((ky(lw(i))/kz(lw(i)))^0.25+(kz(lw(i))/ky(lw(i)))^0.25);
                                          % equivalent radius, ft
    Jw(i) = 2*pi*sqrt(kx(lw(i))*ky(lw(i)))*deltax(lw(i))/(mu*(log(req(i)/rw(i))+S(i)))*6.33e-3;
   end
                                          % Productivity index for BHP well
   if Type(i)==0
        J(lw(i),lw(i)) = Jw(i);
        Wfl(lw(i)) = Jw(i)*BHP(i);
    else
        Wfl(lw(i)) = Rate(i);
    end
end

%***************************BOUNDARY CONDITIONS****************************
Px1 = 5000; % Left Dirichlet boundary pressure, psi
Px2 = 1200; % Right Dirichlet boundary pressure, psi
Py1 = 3000; % Bottom Dirichlet boundary pressure, psi
Py2 = 4000; % Top Dirichlet boundary pressure, psi
wx1 = 0;    % Left boundary BC type (0 - Neumann, 1 - Dirichlet)
wx2 = 0;    % Right boundary BC type
wy1 = 0;    % Bottom boundary BC type
wy2 = 0;    % Top boundary BC type

%***********************T,B,Q,G ASSIGNMENT*********************************
for l = 1:Nx*Ny
    
    if mod(l-1,Nx)~=0              % The grid block is not on left boundary
        kV = [kx(l-1);kx(l)];
        dxV = [deltax(l-1);deltax(l)];
        dyV = [deltay(l-1);deltay(l)];
        hV = [h(l-1);h(l)];
        T(l,l-1) = -trans(kV,Bw,mu,dxV,dyV,hV,3);
        T(l,l) = T(l,l) - T(l,l-1);
     
    elseif wx1 == 1                   % Left face of block has Dirichlet BC
        kV = kx(l);
        dxV = deltax(l);
        dyV = deltay(l);
        hV = h(l);
        T(l,l) = T(l,l) + 2*trans(kV,Bw,mu,dxV,dyV,hV,-3);
        Q(l)= Q(l) + 2*trans(kV,Bw,mu,dxV,dyV,hV,-3)*(Px1-rho*z(l));
    end
   
    if mod(l,Nx)~=0               % The grid block is not on right boundary
        kV = [kx(l);kx(l+1)];
        dxV = [deltax(l); deltax(l+1)];
        dyV = [deltay(l); deltay(l+1)];
        hV = [h(l); h(l+1)];
        T(l,l+1) = -trans(kV,Bw,mu,dxV,dyV,hV,2);
        T(l,l) = T(l,l) - T(l,l+1);
    
    elseif wx2 == 1                  % Right face of block has Dirichlet BC 
        kV = kx(l);
        dxV = deltax(l);
        dyV = deltay(l);
        hV = h(l);
        T(l,l) = T(l,l) + 2*trans(kV,Bw,mu,dxV,dyV,hV,-2);
        Q(l)= Q(l) + 2*trans(kV,Bw,mu,dxV,dyV,hV,-2)*(Px2-rho*z(l));
    end
    
    if l>Nx                      % The grid block is not on bottom boundary
        kV = [ky(l-Nx);ky(l)];
        dxV = [deltax(l-Nx); deltax(l)];
        dyV = [deltay(l-Nx); deltay(l)];
        hV = [h(l-Nx); h(l)];
        T(l,l-Nx) = -trans(kV,Bw,mu,dxV,dyV,hV,4);
        T(l,l) = T(l,l) - T(l,l-Nx);
    
    elseif wy1 == 1                 % Bottom face of block has Dirichlet BC 
        kV = ky(l);
        dxV = deltax(l);
        dyV = deltay(l);
        hV = h(l);
        T(l,l) = T(l,l) + 2*trans(kV,Bw,mu,dxV,dyV,hV,-4);
        Q(l)= Q(l) + 2*trans(kV,Bw,mu,dxV,dyV,hV,-4)*(Py1-rho*z(l));
    end
    
    if l<=Nx*(Ny-1)                 % The grid block is not on top boundary
        kV = [ky(l);ky(l+Nx)];
        dxV = [deltax(l); deltax(l+Nx)];
        dyV = [deltay(l); deltay(l+Nx)];
        hV = [h(l); h(l+Nx)];
        T(l,l+Nx) = -trans(kV,Bw,mu,dxV,dyV,hV,1);
        T(l,l) = T(l,l) - T(l,l+Nx);
    
    elseif wy2 == 1                    % Top face of block has Dirichlet BC 
        kV = ky(l);
        dxV = deltax(l);
        dyV = deltay(l);
        hV = h(l);
        T(l,l) = T(l,l) + 2*trans(kV,Bw,mu,dxV,dyV,hV,-1);
        Q(l)= Q(l) + 2*trans(kV,Bw,mu,dxV,dyV,hV,-1)*(Py2-rho*z(l)); 
    end
    
    B(l,l) = deltax(l)*deltay(l)*h(l)*phi(l)*ct(l)/Bw;
    Q(l) = Q(l) + Wfl(l);
       
end
G = rho*T*z;                                               % Gravity vector
%===============================FLOW SOLVER================================
q = zeros(nt+1,length(lw));
Pwf = zeros(nt,length(lw));
A = T+J+B/dt;
B = sparse(B);
Q = sparse(Q);
A = sparse(A);
Pold = P0;   
tc=[t0;t];
qc=zeros(nt+1,length(lw));

for count=1:nt
    
    t(count) = count*dt+t0;
    tc(count+1) = t(count);
    if t(count)==300                                 % Well Schedule change
        for i=1:length(lw)
            Q(lw(i)) = Q(lw(i)) - Wfl(lw(i));
        end
        Type = [0 0 0 0 0 0 0 0]';
        Rate = [0 0 0 0 0 0 0 0]';
        BHP = [1500 1500 1500 1200 1200 2000 2000 2000]';
        J = sparse(Nx*Ny,Nx*Ny);
        Wfl = sparse(Nx*Ny,1);
        
        for i = 1:length(lw)
            if Type(i)==0
                J(lw(i),lw(i)) = Jw(i);
                Wfl(lw(i)) = Jw(i)*BHP(i);
            else
                Wfl(lw(i)) = Rate(i);
            end
                Q(lw(i)) = Q(lw(i)) + Wfl(lw(i));
        end
        
        A = T+J+B/dt;
        A = sparse(A);
    end
    
    Pnew = A\(B/dt*Pold+Q+G);
    P(:,count) = Pnew;
    
    for i=1:length(lw)
        if Type(i)==0
            q(count+1,i) = Jw(i)*(BHP(i)-P(lw(i),count));
            if count==1
                q(count,i)=q(count+1,i);
            end
            qc(count+1,i)=qc(count,i)+ ((q(count,i)+q(count+1,i))/2)*...
                                             (tc(count+1)-tc(count));
            Pwf(count,i) = BHP(i);
        else
            q(count+1,i) = Rate(i);
            if count==1
                q(count,i)=q(count+1,i);
            end
            qc(count+1,i)=qc(count,i)+ ((q(count,i)+q(count+1,i))/2)*...
                                            (tc(count+1)-tc(count));
            Pwf(count,i) = Rate(i)/Jw(i) + P(lw(i),count);
        end
    end
    
      Pold = Pnew;    
%       prog = [num2str(count/nt*100) '% complete...'];
%       disp(prog)

end

toc
end
