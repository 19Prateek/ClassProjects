clear;
close all;
clc;
%****************************INITIALISATION********************************
%Simulation Properties
Temp = 355.372;                  %K
P=  1.1721e+7;                   %Pa
Tc= [304.2 425.12 617.70 ];       %K
Pc= [ 73.765e5 30.70e5 21.10e5  ];     %Pa
w= [ 0.225 0.2511 0.4898 ];
Rgas= 8.3144598;  %Pa m^3 K^-1 mol^-1
ncomps=3;

%Measured Data
Pmeasured =[7.91 14.80 25.14 35.49 45.83 56.17 62.38 64.51];
xc=[0.765 0.626 0.449 0.330 0.240 0.160 0.116 0.0702];
yc=[0.00466 0.00347 0.00364 0.00493 0.00812 0.0157 0.0306 0.0702];
%Binary Interaction Parameters
Kbinary =[0 0.12 0.1141;
          0.12 0 0
          0.1141 0 0];
aij=Kbinary;
Z=[0,0,0];        %Cubic Roots Initialisation

%Fixed Composition Initialisation
Zfixed= [0.84 0.122 0.038];

%Matrix values initialisation
lnphi=Tc;f=Tc;phi=Tc;phiPure=Tc;
xC1=0.01:0.001:0.999;
delGmixRTMAT=xC1;
DRMAT=xC1;
phiMAT=zeros(size(xC1,2),ncomps);
phiPureMAT=phiMAT;

%***************************PR EOS*****************************************
%Prange=7.91e5:1e5:150e5;%Pa
%for i=1:size(Prange,2)
    %P=Prange(1,i);
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
    flash=0;
    if norm(abs(y-Zfixed))< 1e-9
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
            [iter3,K,vapourFr,liquidFr,V,L] = TwoPFlash(Ki,Pc,Tc,w,P,Kbinary,Temp,Rgas,ncomps,aij,Zfixed,Z);
            nextGuess=0;
            flash=1;
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
        if norm(abs(x-Zfixed))< 1e-9
            disp ('Trivial solution found at composition :' );
            disp (num2str(x));
            disp ('This is second guess for Xi.Stop.');
        elseif x~=Zfixed
            if sum(Xi2)>1
                disp ('Unstable solution found at composition :' );
                disp (num2str(x));
                disp ('Proceed to Two Phase Flash');
                flash=1;
                [iter3,K,vapourFr,liquidFr,V,L] = TwoPFlash(Ki,Pc,Tc,w,P,Kbinary,Temp,Rgas,ncomps,aij,Zfixed,Z);
            elseif sum(Xi2)<1
                disp  ('Solution 2 is stable. Stability established.')
            end
        end
    end
       
    if flash==1
        figure(1)
        xlim([0 1])
        plot(vapourFr(1),P/1e5,'*','Color',[.8;.1;.1])
        hold on
        grid on
        plot(liquidFr(1),P/1e5,'*','Color',[.1;.2;.9])
        %if((P/1e5)>=20 && (P/1e5)<=21.5) || ((P/1e5) >30 &&(P/1e5) <31.5) || ((P/1e5) >40 && (P/1e5)<41.5) || ((P/1e5) >50 && (P/1e5)<51.5)|| ((P/1e5) >60 && (P/1e5)<61.5)
        g=[vapourFr(1),liquidFr(1)];
        h=[P/1e5,P/1e5];
        plot(g,h)
        xlim([0 1])
        disp (['Vapour Fractions =    ' num2str(vapourFr) ]);
        disp (['Liquid Fractions =    ' num2str(liquidFr) ]);
        disp (['V =  ' num2str(V)]);
        disp (['L =  ' num2str(L)]);
        %             disp (['Pressure(bar) =  ' num2str(P/1e5)]);
        disp (['K-Values =  ' num2str(K)]);
        %end
    end
%end

figure (1)
xlabel('Composition CO2')
% ylabel('Pressure(bar)')
% % plot(1-xc,Pmeasured,'X','Color',[.1;.1;.1])
% % plot(1-yc,Pmeasured,'X','Color',[.1;.1;.1])
% legend('Simulation Dew Point Curve', 'Simulation Bubble Point Curve')
