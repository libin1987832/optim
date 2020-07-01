% A=[1,0;-1,0;1,1];
% b=[1,1,3]';
% x0=[10;-10];
% [Q,R]=qr(A'*A);
% AP=pinv(A);
% [x11,f1]=FM(x0,AP,0,A,b)

A=[1,1;-1,-1;-1,0];
b=[10,10,50]';
x0=[-100;-100];
[Q,R]=qr(A);
all=[];
x0GS=x0;
for i=1:10
[xk,r0,rk,fk,fm,fr]=FM(x0,Q,R,A,b);

all=[all [xk;xkGS]];
x0=xk;
x0GS=xkGS;
end
all
[xkGS,r0GS,rkGS,fkGS,fmGS,frGS]=FMGS(x0GS,A,b);