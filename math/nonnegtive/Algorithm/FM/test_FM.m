% A=[1,0;-1,0;1,1];
% b=[1,1,3]';
% x0=[10;-10];
% [Q,R]=qr(A'*A);
% AP=pinv(A);
% [x11,f1]=FM(x0,AP,0,A,b)

% A=[1,1;-1,-1;-1,0];
% b=[10,10,50]';
% x0=[-100;-100];
m=100;n=70;
A=100*rand(m,n)-1;
b=2*rand(m,1)-1;
x0=zeros(n,1);

[Q,R]=qr(A);
all=[];
x0GS=x0;
xk=x0;
for i=1:100


if i>2
    all=[all norm(A*(xk-x0))/norm(A*(x0-x1))];
end
x1=x0;
x0=xk;
    [xk,r0,rk,fk,fm,fr]=FM(x0,Q,R,A,b);

% 
% all=[all [xk;xkGS]];
% x0=xk;
% x0GS=xkGS;
end
all
% all
% [xkGS,r0GS,rkGS,fkGS,fmGS,frGS]=FMGS(x0GS,A,b);