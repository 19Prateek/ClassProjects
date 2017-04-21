function[]= plotting(xpos,ypos,P,Nx,Ny,t,q,Pwf,z)

%=============================PLOTTING=====================================
[X,Y]= meshgrid(xpos,ypos);
earlyt = 10;                                    % Early time
[~,ind_earlyt] = min(abs(t-earlyt));            % Index in t for early time
midt = 100;                                     % Middle time
[~,ind_midt] = min(abs(t-midt));                % Index in t for middle time
latet = 1000;                                   % Late time
[~,ind_latet] = min(abs(t-latet));              % Index in t for late time
ind=[ind_earlyt,ind_midt,ind_latet];
P1 = P(:,ind(1));
P1(z<0.001) = NaN;
P1M = zeros(Ny,Nx);
P2 = P(:,ind(2));
P2(z<0.001) = NaN;
P2M = zeros(Ny,Nx);
P3 = P(:,ind(3));
P3(z<0.001) = NaN;
P3M = zeros(Ny,Nx);

for j=1:Ny
P1M(j,:) = (P1((1+(j-1)*Nx):(Nx+(j-1)*Nx)))';   % early time matrix
P2M(j,:) = (P2((1+(j-1)*Nx):(Nx+(j-1)*Nx)))';   % mid time matrix
P3M(j,:) = (P3((1+(j-1)*Nx):(Nx+(j-1)*Nx)))';   % late time matrix
end

figure(1)
surf(X,Y,P1M,'EdgeColor','none')
xlabel('x-coordinate (ft)')
ylabel('y-coordinate (ft)')
zlabel('Pressure (psi)')
title('Reservoir Pressure Profile (t = 10 days)')
caxis([1600 3500])
colorbar

figure(2)
surf(X,Y,P2M,'EdgeColor','none')
xlabel('x-coordinate (ft)')
ylabel('y-coordinate (ft)')
zlabel('Pressure (psi)')
title('Reservoir Pressure Profile (t = 100 days)')
caxis([1600 3500])
colorbar

figure(3)
surf(X,Y,P3M,'EdgeColor','none')
xlabel('x-coordinate (ft)')
ylabel('y-coordinate (ft)')
zlabel('Pressure (psi)')
title('Reservoir Pressure Profile (t = 1000 days)')
caxis([1600 3500])
colorbar

figure(4)
plot(t,abs(q(2:end,1)),'LineWidth',1.5)
hold on
grid on
grid minor
plot(t,abs(q(2:end,2)),'LineWidth',1.5)
plot(t,abs(q(2:end,3)),'LineWidth',1.5)
plot(t,abs(q(2:end,4)),'LineWidth',1.5)
plot(t,abs(q(2:end,5)),'LineWidth',1.5)
plot(t,abs(q(2:end,6)+q(2:end,7)+q(2:end,8)),'LineWidth',1.5)
hold off
title('Well Rates')
xlabel('Time (days)')
ylabel('Flow Rate (ft3/d)')
legend('Well 1','Well 2','Well 3','Well 4','Well 5','Well 6',0)

figure(5)
plot(t,Pwf(:,1),'LineWidth',1.5)
hold on
grid on
grid minor
plot(t,Pwf(:,2),'LineWidth',1.5)
plot(t,Pwf(:,3),'LineWidth',1.5)
plot(t,Pwf(:,4)','LineWidth',1.5)
plot(t,Pwf(:,5)','LineWidth',1.5)
plot(t,Pwf(:,6),'LineWidth',1.5)
hold off
title('Well BHPs')
xlabel('Time (days)')
ylabel('Bottom Hole Pressure (psi)')
legend('Well 1','Well 2','Well 3','Well 4','Well 5','Well 6',0)

end