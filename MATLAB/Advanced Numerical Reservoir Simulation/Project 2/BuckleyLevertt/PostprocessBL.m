%function [Sw1,Sw2,Sw3] = PostprocessBL(Sw,Nx,dx,dy,Ny,Lx,t)
% Analytical solution is called and plotted to compare with numerical
% simulation
xpos = zeros(Nx,1);                 
ypos = zeros(Ny,1);
for i=1:Nx
    xpos(i)=(i-1)*dx+dx/2;
end
for i=1:Ny
    ypos(i)=(i-1)*dy+dy/2;
end
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^INDEX FINDER FOR TIME^^^^^^^^^^^^^^^^^^
t1 = 56.76;
[~,t1_index] = min(abs(t1-t));
Sw1 = Sw(:,t1_index);
t2 = 113.53;
[~,t2_index] = min(abs(t2-t));
Sw2 = Sw(:,t2_index);
t3 = 188.53;
[~,t3_index] = min(abs(t3-t));
Sw3 = Sw(:,t3_index);
xd = xpos/Lx;

%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ANALYTICAL SOLUTION CALL^^^^^^^^^^^^^^^^^^^^^
[~,~,Sw_an,xd1,xd2,xd3] = fractionflow;
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^POST PROCESSING^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% figure(1)
% plot([1;xd1],[0.2;Sw_an],'-','linewidth',2.5)
% hold on
% plot(xd,Sw1,'-.','linewidth',2.5)
% plot([1;xd2],[0.2;Sw_an],'-','linewidth',2.5)
% plot(xd,Sw2,'-.','linewidth',2.5)
% plot([1;xd3],[0.2;Sw_an],'-','linewidth',2.5)
% plot(xd,Sw3,'-.','linewidth',2.5)
% hold off
% grid on
% title('SATURATION PROFILES AT DIFFERENT DIMENSIONLESS TIMES')
% xlabel('x_{D}')
% ylabel('S_{w}')
% axis([0 1 0 1])
% legend('t_{D} = 0.1 (analytical)','t_{D} = 0.1 (numerical)','t_{D} = 0.2(analytical)','t_{D} = 0.2 (numerical)','t_{D,bt} = 0.33 (analytical)','t_{D,bt} = 0.33 (numerical)',0)

xlsread('Swnumerical.xlsx'); 
a = ans(:,4);
b = ans(:,6);
c = ans(:,9);
d = ans(:,3);


figure(1)
plot([1;xd1],[0.2;Sw_an],'-','linewidth',2.5)
hold on
plot(a,b,'*','linewidth',2.5)
plot([1;xd2],[0.2;Sw_an],'-','linewidth',2.5)
plot(a,c,'*','linewidth',2.5)
plot([1;xd3],[0.2;Sw_an],'-','linewidth',2.5)
plot(a,d,'*','linewidth',2.5)
hold off
grid on
title('SATURATION PROFILES AT DIFFERENT DIMENSIONLESS TIMES')
xlabel('x_{D}')
ylabel('S_{w}')
axis([0 1 0 1])
legend('t_{D} = 0.1 (analytical)','t_{D} = 0.1 (numerical)','t_{D} = 0.2(analytical)','t_{D} = 0.2 (numerical)','t_{D,bt} = 0.33 (analytical)','t_{D,bt} = 0.33 (numerical)',0)




