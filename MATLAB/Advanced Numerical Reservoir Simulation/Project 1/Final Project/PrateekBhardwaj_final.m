% 2-D Reservoir Simulation Project 1
% PGE 392K: Numerical Simulation of Reservoirs
% Prateek Bhardwaj_final

%***********************SIMULATION CONTROL*********************************
[dt,t0,t_fin,nt,rho,zM,kM,phiM,hM,Ny,Nx,mu,Bw] = SimulationControl();
ct = 1e-5*ones(Nx*Ny,1);                    % total compressibility, psi^-1

%***********************INITIALISATION*************************************
[z,kx,phi,h,xpos,ypos,deltax,deltay,dx,dy,T,B,P,t,Q,P0] =Initialise(Nx,...
                                                  Ny,zM,kM,phiM,hM,nt,rho);
 ky = 0.1*kx;                                  % Permeability in y-dir, mD
 kz = kx;                                      % Permeability in y-dir, mD          
%*******************************WELLS**************************************
[J,Wfl,req,Jw,Rate,BHP,Type,lw,S] = Wells(dx,dy,Nx,Ny,kx,ky,kz,deltax,...
                                          deltay,h,mu);
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
A = sparse(T+J+B/dt);
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

        for i = 1:length(lw)
            if Type(i)==0
                J(lw(i),lw(i)) = Jw(i);
                Wfl(lw(i)) = Jw(i)*BHP(i);
            else
                Wfl(lw(i)) = Rate(i);
            end
                Q(lw(i)) = Q(lw(i)) + Wfl(lw(i));
        end
        
        A = (T+J+B/dt);
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

plotting(xpos,ypos,P,Nx,Ny,t,q,Pwf,z)                % Plotting function
