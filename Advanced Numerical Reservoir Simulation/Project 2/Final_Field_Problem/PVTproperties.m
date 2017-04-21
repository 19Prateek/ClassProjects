function [Bofn,muofunc,cofunc] = PVTproperties
% PVT properties, oil pressure dependent viscosity, compressibility and FVF
data=xlsread('PVT.xlsx'); 
Po = data(:,1);
Bo = data(:,3);
muo = data(:,5);
co = data(:,7);
Bofn=fit(Po,Bo,'linearinterp');               
muofunc=fit(Po,muo,'linearinterp');            
cofunc=fit(Po,co,'linearinterp');              
end