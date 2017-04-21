function [x]=history_match
% Read in production history
thist = xlsread('Production History-1.xlsx','A6:A38');
w1hist = xlsread('Production History-1.xlsx','B6:B38');
w2hist = xlsread('Production History-1.xlsx','E6:E38');
w3hist = xlsread('Production History-1.xlsx','H6:H38');
w4hist = xlsread('Production History-1.xlsx','K6:K38');
w5hist = xlsread('Production History-1.xlsx','N6:N38');
w6hist = xlsread('Production History-1.xlsx','Q6:Q38');
w1hist = -w1hist;   % periodic production (bbl)
w2hist = -w2hist;
w3hist = -w3hist;
w4hist = -w4hist;
w5hist = -w5hist;
w6hist = -w6hist;

% Solve for right permeability ratio and compressibility
% %x0 = [0.1 1e-5];
% x0 = [0.2028 5];
% 
% options = optimoptions(@fmincon,'MaxIter',50);
% x = fmincon(@myfun,x0,[],[],[],[],[0.08 0.9],[0.25 5],[],options);
% kratio = x(1);
% comp = x(2);
% disp(kratio)
% disp(comp)
% 
%     function F = myfun(x)
%         
%         [tc,qc] = PrateekBhardwaj_final(x(1),x(2));
%         qc(:,6) = qc(:,6)+qc(:,7)+qc(:,8);
%         qc(:,[7 8]) = [];
%         qc = qc/5.615;
%         index_vector = zeros(length(thist),1);
%         
%         for i=1:length(thist)
%             [~,index_vector(i)] = min(abs(tc-thist(i)));
%         end
% 
%         tc = tc(index_vector);
%         qc = qc(index_vector,:);
%         qper = zeros(size(qc));
%         
% 
%         for i = 2:length(tc)
%             qper(i,:) = qc(i,:)-qc(i-1,:);
%         end
% %         keyboard
%         F = [sum((qper(:,1)-w1hist).^2)+...
%             sum((qper(:,2)-w2hist).^2)+...
%             sum((qper(:,3)-w3hist).^2)+...
%             sum((qper(:,4)-w4hist).^2)+...
%             sum((qper(:,5)-w5hist).^2)+...
%             sum((qper(:,6)-w6hist).^2)];
%         
%         
%     end

[tc,qc] = PrateekBhardwaj_History(0.245,5);
qc(:,6) = qc(:,6)+qc(:,7)+qc(:,8);
qc(:,[7 8]) = [];
qc = qc/5.615;
index_vector = zeros(length(thist),1);

for i=1:length(thist)
    [~,index_vector(i)] = min(abs(tc-thist(i)));
end
tc = tc(index_vector);
qc = qc(index_vector,:);
qper = zeros(size(qc));
for i = 2:length(tc)
    qper(i,:) = qc(i,:)-qc(i-1,:);
end

figure(1)
plot(tc,abs(qper(:,1)),'kd')
hold on
plot(thist,-w1hist,'rd')
plot(tc,abs(qper(:,2)),'ks')
plot(thist,-w2hist,'rs')
plot(tc,abs(qper(:,3)),'kx')
plot(thist,-w3hist,'rx')
plot(tc,abs(qper(:,4)),'kv')
plot(thist,-w4hist,'rv')
plot(tc,abs(qper(:,5)),'k*')
plot(thist,-w5hist,'r*')
plot(tc,abs(qper(:,6)),'ko')
plot(thist,-w6hist,'ro')
hold off
grid on
grid minor
title('History Match')
xlabel('Time (days)')
ylabel('Periodic Production (bbl)')
legend('Simulation Well 1','Field Data Well 1','Simulation Well 2','Field Data Well 2',...
    'Simulation Well 3','Field Data Well 3','Simulation Well 4','Field Data Well 4'...
    ,'Simulation Well 5','Field Data Well 5','Simulation Well 6','Field Data Well 6')

end
