%****************************physical critical properties*****************
%Properties
clc;
clear;
P    = 70.70e6;      %Pa
Temp = 383.15;       %K
Tc   = 871.16;       %K
Pc   = 0.553e+6;     %Pa
w    = 1.4678;
Rgas= 8.3144598;  %Pa m^3 K^-1 mol^-1
Z= [0,0,0];
M= 506.972e-3;         %kg/mol;
%Peng Robinson Equation of State
[a,b,A,B,k,alpha]=PR_para(Pc,Tc,w,P,Temp,Rgas);
%Cordano Cubic Roots for Peng-Robinson
[Z,coeff]=cordano(A,B,Z);
V=Z(1)*Rgas*Temp/P;
dens2=M/V;
disp(['Density at 383.15K and 70.70 Mpa using physical critical properties without volume shift : ' num2str(dens2),' kg/m^3'])


%********************optimized critical properties*************************
%Properties
clear;
P    = 70.70e6;      %Pa
Temp = 383.15;       %K
Tc   = 918.25;       %K
Pc   = 0.916e+6;     %Pa
w    = 1.2184;
Rgas= 8.3144598;  %Pa m^3 K^-1 mol^-1
Z= [0,0,0];
M= 506.972e-3;         %kg/mol;
%Peng Robinson Equation of State
[a,b,A,B,k,alpha]=PR_para(Pc,Tc,w,P,Temp,Rgas);
%Cordano Cubic Roots for Peng-Robinson
[Z,coeff]=cordano(A,B,Z);
V=Z(1)*Rgas*Temp/P;
dens=(P*M)/(Z(1)*Rgas*Temp);
dens2=M/V;
disp(['Density at 383.15K and 70.70 Mpa using optimized critical properties without volume shift : ' num2str(dens2),' kg/m^3'])

%**************physical critical properties with volume shift**************
%Properties
clear;
P    = 70.70e6;      %Pa
Temp = 383.15;       %K
Tc   = 871.16;       %K
Pc   = 0.553e+6;     %Pa
w    = 1.4678;
Rgas= 8.3144598;  %Pa m^3 K^-1 mol^-1
Z= [0,0,0];
M= 506.972e-3;         %kg/mol;
d=2.258;
e=0.1823;
%Peng Robinson Equation of State
[a,b,A,B,k,alpha]=PR_para(Pc,Tc,w,P,Temp,Rgas);
%Cordano Cubic Roots for Peng-Robinson
[Z,coeff]=cordano(A,B,Z);
V=Z(1)*Rgas*Temp/P;
%Volume shift Parameters
s=1-(d/((M*1e3)^e));
c=s*b;
Vmod=V-c;
dens2=M/Vmod;
disp(['Density at 383.15K and 70.70 Mpa using physical critical properties  with volume shift: ' num2str(dens2),' kg/m^3'])