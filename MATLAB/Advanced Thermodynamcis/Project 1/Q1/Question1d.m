clear;
clc;

%Properties
Pi= 1e5;       %Pa
Tempi = 320;   %K
Tc= 370;       %K
Pc= 4.25e+6;   %Pa
w= 0.152;
Rgas= 8.3144598;  %Pa m^3 K^-1 mol^-1
fratio=2;
Z=[0,0,0];
V=0.000091:0.000001:3e-3;
tol =1e-6;
lnphiV=1;
lnphiL=0;
Vliquid=zeros(10,1);
Vvapour=zeros(10,1);
Penv=zeros(10,1);

%**************************************************************************
for i=1:11
    P=Pi;
    lnphiV=1;
    lnphiL=0;
    Temp=Tempi+(i-1)*5;
    while abs(lnphiL-lnphiV) > tol
        %Peng Robinson Equation of State
        [a,b,A,B,k,alpha]=PR_para(Pc,Tc,w,P,Temp,Rgas);
        %Cordano Cubic Roots for Peng-Robinson
        [Z,coeff]=cordano(A,B,Z);
        %Cubic Roots for Compressibility factor
        Zl=min(Z);
        Zv=max(Z);
        %Fugacity Equation
        lnphiV=(Zv-1)-log(Zv-B)-(A/(2*sqrt(2)*B))*log((Zv+(1+sqrt(2))*B)/(Zv+(1-sqrt(2))*B));
        lnphiL=(Zl-1)-log(Zl-B)-(A/(2*sqrt(2)*B))*log((Zl+(1+sqrt(2))*B)/(Zl+(1-sqrt(2))*B));
        phiV=exp(lnphiV);
        phiL=exp(lnphiL);
        P=abs(P*(phiL/phiV));
    end
    disp ([num2str(P) ' Mpa']); 
    Vl=Zl*Rgas*Temp/P;
    Vv=Zv*Rgas*Temp/P;

    Vliquid(i,1)=real(Vl);
    Vvapour(i,1)=real(Vv);
    Penv(i,1)=real(P);

    %PV Isotherms
    V=0.000091:0.000001:Vv+3e-3;
    V=V.';
    aMat=a*ones(size(V));
    bMat=b*ones(size(V));
    TempMat=Temp*ones(size(V));
    PisoT=ones(size(V,1),1);
    PisoT(:,1)=((Rgas*TempMat)./(V-bMat))-(aMat./(V.*(V+bMat)+bMat.*(V-bMat)));
    for j=1:size(V,1)
        if V(j,1)<=Vv && V(j,1)>= Vl
            PisoT(j,1)=P;
        end
    end
    
%Plotting
    figure(1)
    hold on
    semilogx(V,PisoT/1e6,'LineWidth',2)
    ylim([1,7])
    xlim([0,2e-3])
    grid on
    title('PV Isotherms,Critical Point')
    xlabel('Molar Volume(m^3/mol)')
    ylabel('Pressure (MPa)')
    hold on
end
grid minor
legend('320 K','325 K','330 K','335 K','340 K','345 K','350 K','355 K','360 K','365 K','370 K',0)
x1 =  2.2500e-04;
Vliquid(11)=x1;
Vvapour(11)=x1;
y1 = 4.25101;
str1 = '\rightarrow Critical Point';
text(x1,y1,str1,'HorizontalAlignment','left')
hold on
figure(1)
plot(Vliquid,Penv/1e6,'b-')
plot(Vvapour,Penv/1e6,'b-')
