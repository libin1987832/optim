clear
addpath(genpath(pwd));
A=[1,1;-1,-1;-1,0];
b=[1;1;1/2];
x=[-1;-1];
r=b-A*x;
r(r<0)=0;
AA=A'*A;
Ar=A'*r;
AA\Ar
[Q,R]=qr(A);
[x11,f1]=FM(x,Q,R,A,b);