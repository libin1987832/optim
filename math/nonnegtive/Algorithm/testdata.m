% m=2000;
% n=0.1*m;
% A=rand(m,n)*2-ones(m,n);
% b=rand(m,1)*2-ones(m,1);
% x0=zeros(n,1);
% 
% [xk1,fk1,xkArr1]=hybrid1(x0,A,b);
% [xk2,fk2,xkArr2]=hybrid2(x0,A,b);
%[xk3,fk3,xkArr3]=hybrid3(x0,A,b);


% A=[-1,-1;1,1;-1,0;1,0;0,1;0,-1];
% b=[1,1,0.5,-0.5,0.5,-0.5]'
% x0=zeros(2,1);

% x0=-10;
% A=[1;-1];
% b=[100;2];
% [Q,R]=qr(A);
% [xk,r0,rk,fk,fm,fr]=FM(x0,Q,R,A,b)

A=[1,1;-1,-1;-1,0;-6,-3];
b=[1;1;0.5;2];
% x0=[-57/100;47/100];
%  x0=[-3/4;-1/4];
x0=[-3/4;-100];
[Q,R]=qr(A);
[xk,r0,rk,fk,fm,fr]=FM(x0,Q,R,A,b);
r=b-A*xk;
Nk=r;
Nk(Nk>0)=1;
Nk(Nk<0)=0;
for i = 1:20
%     [xk,r0,rk,fk,fm,fr]=FM(x0,Q,R,A,b);
    [xk,r0,rk,z,fk,fr,fu,fz,fd]=FMTD(x0,Q,R,A,b);
    fu/fk
    x0=xk;
    r=b-A*xk;
    Nk2=r;
    Nk2(Nk2>0)=1;
    Nk2(Nk2<0)=0;
    if sum(Nk2==Nk) < 4
        disp(['fd']);
    end
end
r
% diag(ones(4,1))-Q(:,[1,2])*Q(:,[1,2])'*diag([1,1,0,0])
% for i = 1:20
%     [xk,r0,rk,fk,fm,fr]=FM(x0,Q,R,A,b);
% x0=xk; 
% rk
% end
% x0 = [-11/20;11/50];
% A = [1,1;-1,-1;-1,0;-8,-3];
% b = [1;1;1/2;3];
% [xk,fk,y]=ssqr2(x0,A,b)
% b-A*xk
% [xk,fk,y]=ssqr2(xk,A,b)
% b-A*xk
% [xk,fk,y]=ssqr2(xk,A,b)
% b-A*xk
% [xk,fk,y]=ssqr2(xk,A,b)
% b-A*xk



% A = [-1,-1;1,1;-1,0];
% b = [1,1,0.5]';
% x0 = [-10,-10]';
% 
% ATA = A'*A;
% [v,d]=eig(ATA);
% lambal = max(max(d));
% for i = 1:30
%     [xk,r0,rk,fk,fm,fr] = BFM(x0,lambal,A,b);
%     x0 = xk
% end
% [xk1,fk1,xkArr1]=hybrid1(x0,A,b);
% [xk2,fk2,xkArr2]=hybrid2(x0,A,b);
% [xk3,fk3,xkArr3]=hybrid3(x0,A,b);