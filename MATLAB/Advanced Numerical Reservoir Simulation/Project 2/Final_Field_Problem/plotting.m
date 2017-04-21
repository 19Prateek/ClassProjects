function[]= plotting(dx,dy,Po,Po0,Sw,Sw0,Nx,Ny,t,qo,qw,Pwf,z)
xpos = zeros(Nx,1);                
ypos = zeros(Ny,1);
for i=1:Nx
    xpos(i)=(i-1)*dx+dx/2;
end

for i=1:Ny
    ypos(i)=(i-1)*dy+dy/2;
end

q= -(qw+qo)/5.615; % Postive rates are production and negative is injection

[X,Y]= meshgrid(xpos,ypos);
earlyt = 500;                                   
[~,ind_earlyt] = min(abs(t-earlyt));           
midt = 1826;                                    
[~,ind_midt] = min(abs(t-midt));                
latet = 3987;                                   
[~,ind_latet] = min(abs(t-latet));             
ind=[ind_earlyt,ind_midt,ind_latet];
P0 = Po0;
P0(z<0.001) = NaN;
P0M = zeros(Ny,Nx);
S0 = Sw0;
S0(z<0.001) = NaN;
S0M = zeros(Ny,Nx);
P1 = Po(:,ind(1));            
P1(z<0.001) = NaN;
P1M = zeros(Ny,Nx);
S1 = Sw(:,ind(1));
S1(z<0.001) = NaN;
S1M = zeros(Ny,Nx);
P2 = Po(:,ind(2));
P2(z<0.001) = NaN;
P2M = zeros(Ny,Nx);
S2 = Sw(:,ind(2));
S2(z<0.001) = NaN;
S2M = zeros(Ny,Nx);
P3 = Po(:,ind(3));
P3(z<0.001) = NaN;
P3M = zeros(Ny,Nx);
S3 = Sw(:,ind(3));
S3(z<0.001) = NaN;
S3M = zeros(Ny,Nx);

for j=1:Ny
P0M(j,:) = (P0((1+(j-1)*Nx):(Nx+(j-1)*Nx)))';   % Pressure at t=0 days
P1M(j,:) = (P1((1+(j-1)*Nx):(Nx+(j-1)*Nx)))';   % Pressure at t=500 days
P2M(j,:) = (P2((1+(j-1)*Nx):(Nx+(j-1)*Nx)))';   % Pressure at t=1826 days
P3M(j,:) = (P3((1+(j-1)*Nx):(Nx+(j-1)*Nx)))';   % Pressure at t=3897 days
S0M(j,:) = (S0((1+(j-1)*Nx):(Nx+(j-1)*Nx)))';   % Saturation at t=0 days
S1M(j,:) = (S1((1+(j-1)*Nx):(Nx+(j-1)*Nx)))';   % Saturation at t=500 days
S2M(j,:) = (S2((1+(j-1)*Nx):(Nx+(j-1)*Nx)))';   % Saturation at t=1826 days
S3M(j,:) = (S3((1+(j-1)*Nx):(Nx+(j-1)*Nx)))';   % Saturation at t=3897days
end

 % Initial time Pressure profile(t=0 days)
figure(1)
surf(X,Y,P0M,'EdgeColor','none')
xlabel('x-coordinate (ft)')
ylabel('y-coordinate (ft)')
zlabel('Pressure (psi)')
title('Reservoir Oil Pressure Profile (t = 0 days)')
colorbar

 % Early time Pressure profile(t=500 days)
figure(2)
surf(X,Y,P1M,'EdgeColor','none')
xlabel('x-coordinate (ft)')
ylabel('y-coordinate (ft)')
zlabel('Pressure (psi)')
title('Reservoir Oil Pressure Profile (t = 500 days)')
colorbar

 % Middle time Pressure profile(t=1826 days)
figure(3)
surf(X,Y,P2M,'EdgeColor','none')
xlabel('x-coordinate (ft)')
ylabel('y-coordinate (ft)')
zlabel('Pressure (psi)')
title('Reservoir Oil Pressure Profile (t = 1826 days)')
%caxis([1600 3500])
colorbar

% Late Time Pressure profile(t=3897 days)
figure(4)
surf(X,Y,P3M,'EdgeColor','none')
xlabel('x-coordinate (ft)')
ylabel('y-coordinate (ft)')
zlabel('Pressure (psi)')
title('Reservoir Oil Pressure Profile (t = 3987 days)')
colorbar

%Initial Time Saturation profile(t=0 days)
figure(5)
surf(X,Y,S0M,'EdgeColor','none')
xlabel('x-coordinate (ft)')
ylabel('y-coordinate (ft)')
zlabel('Saturation')
title('Reservoir Water Saturation Profile (t = 0 days)')
colorbar

%Early Time Saturation profile(t=500 days)
figure(6)
surf(X,Y,S1M,'EdgeColor','none')
xlabel('x-coordinate (ft)')
ylabel('y-coordinate (ft)')
zlabel('Saturation')
title('Reservoir Water Saturation Profile (t = 500 days)')
colorbar

% Middle Time Saturation profile(t=1826 days)
figure(7)
surf(X,Y,S2M,'EdgeColor','none')
xlabel('x-coordinate (ft)')
ylabel('y-coordinate (ft)')
zlabel('Saturation')
title('Reservoir Water Saturation Profile (t = 1826  days)')
colorbar

% Late Time Saturation profile(t=3897 days)
figure(8)
surf(X,Y,S3M,'EdgeColor','none')
xlabel('x-coordinate (ft)')
ylabel('y-coordinate (ft)')
zlabel('Saturation ')
title('Reservoir Water Saturation Profile (t = 3987 days)')
colorbar

figure(9)
plot(t,q(2:end,1),'LineWidth',1.5)
hold on
grid minor
plot(t,q(2:end,2),'LineWidth',1.5)
plot(t,q(2:end,3)+q(2:end,4)+q(2:end,5),'LineWidth',1.5)
plot(t,q(2:end,6)+q(2:end,7)+q(2:end,8),'LineWidth',1.5)
plot(t,q(2:end,9),'LineWidth',1.5)
plot(t,q(2:end,10),'LineWidth',1.5)
hold off
title('Well Rates')
xlabel('Time (days)')
ylabel('Flow Rate (bbl/d)')
legend('Well 1','Well 2','Well 3','Well 4','Well 5','Well 6',0)

Pwf(1826:1830,9)=Pwf(1831,9);
Pwf(1826:1830,10)=Pwf(1831,10);

figure(10)
plot(t,Pwf(:,1),'LineWidth',1.5)
hold on
grid minor
plot(t,Pwf(:,2),'LineWidth',1.5)
plot(t,(Pwf(:,3)+Pwf(:,4)+Pwf(:,5))/3,'LineWidth',1.5)
plot(t,(Pwf(:,6)+Pwf(:,7)+Pwf(:,8))/3,'LineWidth',1.5)
plot(t,Pwf(:,9)','LineWidth',1.5)
plot(t,Pwf(:,10),'LineWidth',1.5)
hold off
title('Well BHPs')
xlabel('Time (days)')
ylabel('Bottom Hole Pressure (psi)')
legend('Well 1','Well 2','Well 3','Well 4','Well 5','Well 6',0)

end