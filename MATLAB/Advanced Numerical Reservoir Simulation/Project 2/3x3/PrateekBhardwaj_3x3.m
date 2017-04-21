clear
clc
[reservoir,numerical,well] = Initialise();
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^INITIALISATION^^^^^^^^^^^^^^^^^^^^^^^^^^^
dt = numerical.dt;                               
t0 = numerical.t0;                               
t_fin = numerical.t_fin;                         
Nx = reservoir.Nx;
Ny = reservoir.Ny;
muw = reservoir.muw;                               
muo = reservoir.muo;
Bw = reservoir.Bw;                                 
Bo = reservoir.Bo;
cw = reservoir.cw;                                 
co = reservoir.co;                                
cr = reservoir.cr*ones(Nx*Ny,1);                   
SwM = reservoir.Sw0*ones(Ny,Nx);
SoM = 1-SwM;
Lx = reservoir.Lx;                              
Ly = reservoir.Ly;                               
kM = reservoir.kx*ones(Ny,Nx);
phiM = reservoir.phi*ones(Ny,Nx);
hM = reservoir.h*ones(Ny,Nx);
zM = reservoir.z*ones(Ny,Nx);
nt = round((t_fin-t0)/dt);             
Sw0 = zeros(Nx*Ny,1);                  
So0 = 1-Sw0;                         
kx = zeros(Nx*Ny,1);                    
phi = zeros(Nx*Ny,1);                   
h = zeros(Nx*Ny,1);                     
for j=1:Ny
%    z((1+(j-1)*Nx):(Nx+(j-1)*Nx))=(zM(j,:))';
    Sw0((1+(j-1)*Nx):(Nx+(j-1)*Nx))=(SwM(j,:))';
    So0((1+(j-1)*Nx):(Nx+(j-1)*Nx))=(SoM(j,:))';
    kx((1+(j-1)*Nx):(Nx+(j-1)*Nx))=(kM(j,:))';
    phi((1+(j-1)*Nx):(Nx+(j-1)*Nx))=(phiM(j,:))';
    h((1+(j-1)*Nx):(Nx+(j-1)*Nx))=(hM(j,:))';
end
% z=-z;                           
ky = kx;                         
kz = kx;                            
dx = Lx/Nx;                      
dy = Ly/Ny;                        
deltax = dx*ones(Nx*Ny,1);         
deltay = dy*ones(Nx*Ny,1);         
reservoir.kx = kx;
reservoir.ky = ky;
reservoir.kz = kz;
reservoir.h = h;
reservoir.phi = phi;
reservoir.deltay = deltay;
reservoir.deltax = deltax;
reservoir.dx = dx;
reservoir.dy = dy;
numerical.nt = nt;
% zmax = max(z);
P0 = 1000*ones(Nx*Ny,1);                

%^^^^^^^^^^^^^^^^^^^^^^^^^^PREALLOCATE^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
P = zeros(Nx*Ny,nt);                   % P matrix for recording pressures
Sw = zeros(Nx*Ny,nt);                  % Sw matrix
So = zeros(Nx*Ny,nt);
t = zeros(nt,1);                       % t vector for recording times
Qw = sparse(Nx*Ny,1);                  % Qw matrix used in several functions
Qo = Qw;
well.Wflw = zeros(length(well.S),1);   % Well flow rate (water) vector to be used
well.Wflo = zeros(length(well.S),1);   % Well flow rate (oil) vector to be used
well.nS = 1;                           % well schedule hasn't been assigned

%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^FLOW SOLVER^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
q = zeros(nt+1,length(well.S));
Pwf = zeros(nt,length(well.S));
tc=[t0;t];
qc=zeros(nt+1,length(well.S));
Pold = P0;
Swold = Sw0;
SC = 0;           

for count=1:nt
    
    t(count) = count*dt+t0;
    tc(count+1) = t(count);  
    [T,Tw,To,D,reservoir] = TDG(numerical,reservoir,Pold,Swold);
    [well,Qw,Qo,Q,J,Jwater] = Wells(reservoir,well,Qw,Qo,Swold,SC);
    
    A = T+J+D;
    A = sparse(A);
    if well.nS<well.nSC+1
        if t(count)==well.tSC(well.nS)            
            SC = 1;
            [well,Qw,Qo,Q,J,Jwater] = Wells(reservoir,well,Qw,Qo,Swold,SC); 
            A = T+J+D;
            A = sparse(A);
        end
    end
    Pnew = A\(D*Pold+Q);
    Swnew = Swold+(reservoir.d12)^(-1)*(-reservoir.d11*(Pnew-Pold)+(Qw-Jwater*Pnew)-Tw*(Pnew));
    P(:,count) = Pnew;
    Sw(:,count) = Swnew;
    Pold = Pnew;
    Swold = Swnew;
    prog = [num2str(count/nt*100) '% complete...'];
    disp(prog)
    
end

xpos = zeros(Nx,1);                   
ypos = zeros(Ny,1);
xpos(1) = 0;
ypos(1) = 0;
for i=1:Nx
    xpos(i)=(i-1)*dx+dx/2;
end

for i=1:Ny
    ypos(i)=(i-1)*dy+dy/2;
end