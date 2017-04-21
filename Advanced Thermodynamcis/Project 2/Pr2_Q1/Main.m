clear;
close all;
clc;
%****************************INITIALISATION********************************
%Simulation Properties
Temp = 621.8 ;                  %K
P=40.53e5;                     %Pa
Tc= [507.60 871.16 ];       %K
Pc= [30.25e5 5.53e5 ];     %Pa
w= [0.3010 1.4678];
Rgas= 8.3144598;  %Pa m^3 K^-1 mol^-1
ncomps=2;

%Binary Interaction Parameters
Kbinary =[0 0;
    0 0];
aij=Kbinary;
Z=[0,0,0];        %Cubic Roots Initialisation

%Fixed Composition Initialisation
Zfixed= [0.9 0.1];

%Matrix values initialisation
lnphi=Tc;f=Tc;phi=Tc;phiPure=Tc;
xC1=0.01:0.001:0.999;
delGmixRTMAT=xC1;
DRMAT=xC1;
phiMAT=zeros(size(xC1,2),ncomps);
phiPureMAT=phiMAT;

%***************************PR EOS*****************************************
%Peng Robinson EOS parameters
[a,b,A,B,~]=PR_para(Pc,Tc,w,P,Temp,Rgas,ncomps);

%********************************Fixed Composition analysis****************
[phiZfixed,Zfixed] = phiFixed(Zfixed,b,a,B,Kbinary,aij,P,Rgas,Temp,Tc,ncomps,Z);

%Wilson's Correlation for intial estimate
Ki=(Pc/P).*exp(5.373.*(1+w).*(1-(Tc/Temp)));

%*****************************Vapour Like Estimate*************************
Xi=Zfixed.*Ki;
[y,Xi,iter1] = stationaryPoint(Xi,Zfixed,b,a,B,aij,Kbinary,Tc,Rgas,Temp,P,ncomps,Z,phiZfixed);

disp ('Zfixed :' );
disp (num2str(Zfixed));

if norm(abs(y-Zfixed))< 1e-8
    disp ('Trivial solution found at composition :' );
    disp (num2str(y));
    disp ('Moving to second guess for Xi.');
    Xi=Zfixed./Ki;
    [x,Xi2,iter2] = stationaryPoint(Xi,Zfixed,b,a,B,aij,Kbinary,Tc,Rgas,Temp,P,ncomps,Z,phiZfixed);
    nextGuess=1;
elseif y~=Zfixed
    if sum(Xi)>1
        disp ('Unstable solution found at composition :' );
        disp (num2str(y));
        disp ('Proceed to Two Phase Flash');
        [iter3,K] = TwoPFlash(Ki,Pc,Tc,w,P,Kbinary,Temp,Rgas,ncomps,aij,Zfixed,Z);
        nextGuess=0;
    elseif sum(Xi)<1
        disp ('Stable solution found at composition :' );
        disp (num2str(y));
        disp ('Checking at second guess for Xi.');
        Xi=Zfixed./Ki;
        nextGuess=1;
        [x,Xi2,iter2] = stationaryPoint(Xi,Zfixed,b,a,B,aij,Kbinary,Tc,Rgas,Temp,P,ncomps,Z,phiZfixed);
    end
end

if nextGuess==1
    if norm(abs(x-Zfixed))< 1e-8
        disp ('Trivial solution found at composition :' );
        disp (num2str(x));
        disp ('This is second guess for Xi.Stop.');
    elseif x~=Zfixed
        if sum(Xi2)>1
            disp ('Unstable solution found at composition :' );
            disp (num2str(x));
            disp ('Proceed to Two Phase Flash');
            [iter3,K] = TwoPFlash(Ki,Pc,Tc,w,P,Kbinary,Temp,Rgas,ncomps,aij,Zfixed,Z);
        elseif sum(Xi2)<1
            disp  ('Solution 2 is stable. Stability established.')
        end
    end
end





