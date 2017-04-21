function [z,kx,phi,h,xpos,ypos,deltax,deltay,dx,dy,T,B,P,t,Q,P0] =Initialise(Nx,Ny,zM,kM,phiM,hM,nt,rho)

%******************************MATRIX TO VECTORS***************************
z = zeros(Nx*Ny,1);                     % Depth , ft
kx = zeros(Nx*Ny,1);                    % Permeability in x-dir, mD
phi = zeros(Nx*Ny,1);                   % Porostiy,fraction
h = zeros(Nx*Ny,1);                     % Thickness, ft

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

xpos(1) = 0;
ypos(1) = 0;
for i=1:Nx
    xpos(i)=(i-1)*dx+dx/2;
end

for i=1:Ny
    ypos(i)=(i-1)*dy+dy/2;
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

end

