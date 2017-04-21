function [J,Wfl,req,Jw,Rate,BHP,Type,lw,S] = Wells(dx,dy,Nx,Ny,kx,ky,kz,deltax,deltay,h,mu)

Location = [4184 1569; 
            4968.5 2510.4; 
            3294.9 2928.8; 
            2562.7 4393.2; 
            1307.5 2824.2;
            890 895; 
            965 895; 
            815 895];

iwell = ceil(Location(:,1)/dx);  
jwell = ceil(Location(:,2)/dy);
lw = iwell + (jwell - 1) * Nx;             % l indices from x, y coordinate          


Type = [ 1 1 1 0 0 1 1 1]';                % Type 1= constant rate, Type 0= constant BHP
Orientation = [ 0 0 0 0 0 1 1 1]';         %  1=Horizontal Well,  0=Vertical Well
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
                                         
   if Type(i)==0
        J(lw(i),lw(i)) = Jw(i);
        Wfl(lw(i)) = Jw(i)*BHP(i);
    else
        Wfl(lw(i)) = Rate(i);
    end
end


end

