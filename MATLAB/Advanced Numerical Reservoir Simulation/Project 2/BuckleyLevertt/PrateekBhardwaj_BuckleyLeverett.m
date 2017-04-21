clear
clc
[reservoir,numerical,well] = InitialiseBL();

%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^INITIALISATION^^^^^^^^^^^^^^^^^^^^^^^^^^^
Lx = reservoir.Lx;                      
Ly = reservoir.Ly;  
Nx = reservoir.Nx;
Ny = reservoir.Ny;
kM = reservoir.kx*ones(Ny,Nx);
phiM = reservoir.phi*ones(Ny,Nx);
hM = reservoir.h*ones(Ny,Nx);
zM = reservoir.z*ones(Ny,Nx);
dt = numerical.dt;                        
t0 = numerical.t0;                        
t_fin = numerical.t_fin;                  
muw = reservoir.muw;                     
muo = reservoir.muo;
Bw = reservoir.Bw;                       
Bo = reservoir.Bo;
cw = reservoir.cw;                        
co = reservoir.co;                       
cr = reservoir.cr*ones(Nx*Ny,1);          
SwM = reservoir.Sw0*ones(Ny,Nx);
SoM = 1-SwM;
nt = round((t_fin-t0)/dt);                
Sw0 = zeros(Nx*Ny,1);                    
Pc0 = zeros(Nx*Ny,1);
kx = zeros(Nx*Ny,1);                      
phi = zeros(Nx*Ny,1);                     
h = zeros(Nx*Ny,1);                       
for j=1:Ny
    z((1+(j-1)*Nx):(Nx+(j-1)*Nx))=(zM(j,:))';
    Sw0((1+(j-1)*Nx):(Nx+(j-1)*Nx))=(SwM(j,:))';
    kx((1+(j-1)*Nx):(Nx+(j-1)*Nx))=(kM(j,:))';
    phi((1+(j-1)*Nx):(Nx+(j-1)*Nx))=(phiM(j,:))';
    h((1+(j-1)*Nx):(Nx+(j-1)*Nx))=(hM(j,:))';
end
kx(1) = 1000;
kx(end) = 1000;
ky = kx;                            
kz = kx;                            
dx = Lx/Nx;                        
dy = Ly/Ny;                        
deltax = dx*ones(Nx*Ny,1);          
deltay = dy*ones(Nx*Ny,1);          
reservoir.kx = kx;
reservoir.ky = ky;
reservoir.kz = kz;
reservoir.z = z;
reservoir.h = h;
reservoir.phi = phi;
reservoir.deltay = deltay;
reservoir.deltax = deltax;
reservoir.dx = dx;
reservoir.dy = dy;
numerical.nt = nt;
P0 = 1000*ones(Nx*Ny,1);                 

%^^^^^^^^^^^^^^^^^^^^^^^^^^PREALLOCATE^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Po = zeros(Nx*Ny,nt);                  
Pc = zeros(Nx*Ny,nt);
Sw = zeros(Nx*Ny,nt);                 
So = zeros(Nx*Ny,nt);
t = zeros(nt,1);                      
Qw = sparse(Nx*Ny,1);                 
Qo = Qw;
well.Wflw = zeros(length(well.S),1);   
well.Wflo = zeros(length(well.S),1);   
well.nS = 1;                         

%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^FLOW SOLVER^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
q = zeros(nt+1,length(well.S));
Pwf = zeros(nt,length(well.S));
tc=[t0;t];
qc=zeros(nt+1,length(well.S));
Pold = P0;
Swold = Sw0;
Pcold = zeros(Nx*Ny,1);
SC = 0;                        
for count=1:nt
    t(count) = count*dt+t0;
    tc(count+1) = t(count);
    [T,Tw,To,D,G,reservoir] = TDGBL(numerical,reservoir,Pold,Swold,Pcold);
    [well,Qw,Qo,Q,J] = WellsBL(reservoir,well,Qw,Qo,Swold,Pcold,SC);
    A = T+J+D;
    A = sparse(A);
    if well.nS<well.nSC+1
        if t(count)==well.tSC(well.nS)            
            SC = 1;
            [well,Qw,Qo,Q,J] = WellsBL(reservoir,well,Qw,Qo,Swold,Pcold,SC); 
            A = T+J+D;
            A = sparse(A);
        end 
    end
    Pnew = A\(D*Pold+Q+G);
    Swnew = Swold+(reservoir.d12)^(-1)*(-reservoir.d11*(Pnew-Pold)+(Qw-well.Jwater*Pnew)-Tw*(Pnew-0));
    Pcnew = zeros(Nx*Ny,1);
    Po(:,count) = Pnew;
    Sw(:,count) = Swnew;
    Pold = Pnew;
    Swold = Swnew;
    Pcold = Pcnew;
    prog = [num2str(count/nt*100) '% complete...'];
    disp(prog)
    
end

PostprocessBL(Sw,Nx,dx,dy,Ny,Lx,t)


