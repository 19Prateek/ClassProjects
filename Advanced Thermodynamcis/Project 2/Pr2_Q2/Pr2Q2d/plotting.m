%Plottin
close all
figure(3)
plot(xC1,delGmixRTMAT,'.')
%ylim([-50,30])
xlim([0,1])
grid on
grid minor
% title('PV Isotherm, T=383.15K deg C')
xlabel('x1(fraction)')
ylabel('delGmixRT')
hold on
plot(xC1,delGmixRTMAT1,'.')
legend('BIP=0','BIP=0.099',0)

