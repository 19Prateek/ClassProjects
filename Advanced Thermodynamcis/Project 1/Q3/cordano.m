function [Z,coeff]=cordano(A,B,Z)
a2= -(1-B);
a1=  (A-3*B^2-2*B);
a0= -(A*B-B^2-B^3);
coeff=[1,a2,a1,a0];
Q=(3*a1-a2^2)/9;
R =(9*a2*a1-27*a0-2*a2^3)/54;
D=Q^3+R^2;
if(D<0)
    theta=acos(R/sqrt(-Q^3));
    Z(1)=2*sqrt(-Q)*cos(theta/3)-a2/3;
    Z(2)=2*sqrt(-Q)*cos((theta+2*pi())/3)-a2/3;
    Z(3)=2*sqrt(-Q)*cos((theta+4*pi())/3)-a2/3;
elseif(D>0)
    S=nthroot(R+sqrt(D),3);
    T=nthroot(R-sqrt(D),3);
    Z(1)=-a2/3 +(S+T);
    Z(2)=-0.5*(S+T)-(1/3)*a2 + 0.5i*sqrt(3)*(S-T);
    Z(3)=-0.5*(S+T)-(1/3)*a2 - 0.5i*sqrt(3)*(S-T);
end
Z=[Z(1),Z(2),Z(3)]; 