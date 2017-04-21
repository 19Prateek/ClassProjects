% Numerical Reservoir Simulation Project 2
% Contributors : Prateek Bhardwaj, Sai Uppati
clear
clc
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^READ FROM USER INPUTS^^^^^^^^^^^^^^^^^^^^^
% All these inputs have been described and are being read from the input file
[reservoir,numerical,well] = UserInput();

dt = numerical.dt;
t0 = numerical.t0;
t_fin = numerical.t_fin;
Nx = reservoir.Nx;
Ny = reservoir.Ny;
Lx = reservoir.Lx;
Ly = reservoir.Ly;
muw = reservoir.muw;
Bw = reservoir.Bw;
cw = reservoir.cw;
cr = reservoir.cr;
rhow = reservoir.rhow;
rhoo = reservoir.rhoo;
kM = reservoir.kM;
phiM = reservoir.phiM;
hM = reservoir.hM;
zM = reservoir.zM;
nt = round((t_fin-t0)/dt);

%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^MATRIX TO VECTOR^^^^^^^^^^^^^^^^^^^^^^^^^^^^
z = zeros(Nx*Ny,1);
kx = zeros(Nx*Ny,1);
phi = zeros(Nx*Ny,1);
h = zeros(Nx*Ny,1);
for j=1:Ny
    z((1+(j-1)*Nx):(Nx+(j-1)*Nx))=(zM(j,:))';
    %     Sw0((1+(j-1)*Nx):(Nx+(j-1)*Nx))=(SwM(j,:))';    % Buckley Leverette
    kx((1+(j-1)*Nx):(Nx+(j-1)*Nx))=(kM(j,:))';
    phi((1+(j-1)*Nx):(Nx+(j-1)*Nx))=(phiM(j,:))';
    h((1+(j-1)*Nx):(Nx+(j-1)*Nx))=(hM(j,:))';
end
z = -z;
ky = 0.15*kx;
kz = kx;
dx = Lx/Nx;
dy = Ly/Ny;
deltax = dx*ones(Nx*Ny,1);
deltay = dy*ones(Nx*Ny,1);
numerical.nt = nt;
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

%^^^^^^^^^^^^^^^SATURATION INITIALISATION FROM CAPILLARY PRESSURE^^^^^^^^^^
Powoc = reservoir.Powoc;
Pe = 3.5;
Pwwoc = Powoc - Pe;
zwoc = 7474.45;
Po0 = Powoc + rhoo*(z-zwoc);
Pw0 = Pwwoc + rhow*(z-zwoc);
Pc0 = Po0 - Pw0;
Sw0 = pcap(Pc0);
Sw0(z>=zwoc) = 1;

%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^PVT PROPERTIES^^^^^^^^^^^^^^^^^^^^^^^^^^
[~,muofunc,cofunc] = PVTproperties;

%^^^^^^^^^^^^^^^^^^^^^^^^MATRIX & VECTOR PREALLOCATION^^^^^^^^^^^^^^^^^^^^^
Po = zeros(Nx*Ny,nt);
Pc = zeros(Nx*Ny,nt);
Sw = zeros(Nx*Ny,nt);
t = zeros(nt,1);
Qw = sparse(Nx*Ny,1);
Qo = Qw;
well.Wflw = zeros(length(well.S),1);
well.Wflo = zeros(length(well.S),1);
well.nS = 1;
qo = zeros(nt+1,length(well.S));
qw = zeros(nt+1,length(well.S));
Pwf = zeros(nt,length(well.S));
tc=[t0;t];
qoc=zeros(nt+1,length(well.S));
qwc=zeros(nt+1,length(well.S));

Pold = Po0;
Swold = Sw0;
Pcold = Pc0;
Pcold(Pcold>1000) = 1000;
reservoir.BoV = ones(6000,1);

% Time Loop
for count=1:nt
    t(count) = count*dt+t0;
    tc(count+1) = t(count);
    reservoir.muoV = muofunc(Pold);
    reservoir.coV = cofunc(Pold);
    
    %^^^^^^^^^^^^^^^^^^^^^^^TDGQ ASSIGNMENT^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    [T,Tw,To,D,G,reservoir] = TDG(numerical,reservoir,Pold,Swold,Pcold);
    [well,Qw,Qo,Q,J] = Wells(reservoir,well,Qw,Qo,Swold,Pcold);
    
    if well.nS<well.nSC+1
        if t(count)==well.tSC(well.nS)
            well.nS = well.nS + 1;
            [well,Qw,Qo,Q,J] = Wells(reservoir,well,Qw,Qo,Swold,Pcold);
        end
    end
    
    %^^^^^^^^^^^^^^^^^^^^^^^MULTIPHASE FLOW EQUATION LINEAR SOLVER^^^^^^^^^^
    A = T+J+D;
    A = sparse(A);
    Pnew = A\(D*Pold+Q+G);
    Swnew = Swold+(reservoir.d12inv)*(-reservoir.d11*(Pnew-Pold)+ ...
        (Qw-well.Jwater*(Pnew-Pcold))-Tw*(Pnew-Pcold-rhow*z));
    
    %^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^RATES & BHP^^^^^^^^^^^^^^^^^^^^^^^^^^
    Swnew(Swnew<=0.2) = 0.2;
    Swnew(Swnew>=1)= 1;
    for i = 1:length(well.lw)
        Swnew(Swnew(well.lw(i))>0.8)=0.8;
    end
    Pcnew = capillary(Swnew);
    Pcnew(Pcnew>1000) = 1000;
    Po(:,count) = Pnew;
    Sw(:,count) = Swnew;
    Pc(:,count) = Pcnew;
    
    for i=1:length(well.lw)
        if well.Type(i,well.nS)==0
            qo(count+1,i) = well.Jo(i)*(well.BHP(i,well.nS)-Po(well.lw(i),count));
            qw(count+1,i) = well.Jw(i)*(well.BHP(i,well.nS)-(Po(well.lw(i),count)-Pcold(well.lw(i))));
            if count==1
                qo(count,i)=qo(count+1,i);
                qw(count,i)=qw(count+1,i);
            end
            qoc(count+1,i) = qoc(count,i)+ ((qo(count,i)+qo(count+1,i))/2)*...
                (tc(count+1)-tc(count));
            qwc(count+1,i) = qwc(count,i)+ ((qw(count,i)+qw(count+1,i))/2)*...
                (tc(count+1)-tc(count));
            Pwf(count,i) = well.BHP(i,well.nS);
        else
            if well.Rate(i,well.nS) < 0
                qo(count+1,i) = well.fo(i)*well.Rate(i,well.nS);
                qw(count+1,i) = well.fw(i)*well.Rate(i,well.nS);
                if count==1
                    qo(count,i)=qo(count+1,i);
                    qw(count,i)=qw(count+1,i);
                end
                qoc(count+1,i)=qoc(count,i)+ ((qo(count,i)+qo(count+1,i))/2)*(tc(count+1)-tc(count));
                qwc(count+1,i)=qwc(count,i)+ ((qw(count,i)+qw(count+1,i))/2)*(tc(count+1)-tc(count));
                Pwf(count,i) = well.fo(i)*well.Rate(i,well.nS)/well.Jo(i)+  (Po(well.lw(i),count));
            else
                qo(count+1,i) = 0;
                qw(count+1,i) = well.Rate(i,well.nS);
                if count==1
                    qo(count,i)=qo(count+1,i);
                    qw(count,i)=qw(count+1,i);
                end
                qoc(count+1,i)=qoc(count,i)+ ((qo(count,i)+qo(count+1,i))/2)*(tc(count+1)-tc(count));
                qwc(count+1,i)=qwc(count,i)+ ((qw(count,i)+qw(count+1,i))/2)*(tc(count+1)-tc(count));
                Pwf(count,i) = (well.Rate(i,well.nS))/(well.Jw(i)) + (Po(well.lw(i),count)-Pcold(well.lw(i)));
            end
        end
    end
    Pold =   Pnew;
    Swold =  Swnew;
    Pcold =  Pcnew;
    prog = [num2str(count/nt*100) '% complete...'];
    disp(prog)
    
end

%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^PLOTS/POST PROCESSING^^^^^^^^^^^^^^^^^^^^^^
plotting(dx,dy,Po,Po0,Sw,Sw0,Nx,Ny,t,qo,qw,Pwf,z)
