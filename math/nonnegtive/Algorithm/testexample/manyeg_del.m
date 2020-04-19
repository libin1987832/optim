% A=[1 ;-1 ;0 ];
% b=[2;3;1];
% x0=[1];
% [q,r]=qr(A);
% [x1,f1]=FM(x0,q,r,A,b)
% dFM(A,b,x1)
A=[1,1;-1,-1;-1,0;-8,-3]
b=[1;1;1/2;3];
x0=[-0.50000000бо1;0.3333333];
[q,r]=qr(A);
q3=q(3,[1,2])*q(:,[1,2])';
y0=b-A*x0;
yp=y0;
yp(y0<0)=0;
y1=y0(3)-q3*yp;
[xk,r0,rk,fk,fm,fr]=FM(x0,q,r,A,b);
y1x=b-A*xk;
y1
y1x(3)
% A=[1,1;-1,-1;-66,-60;-1,0];
% b=[1;1;33;1/2];
% x0=[-1;1/2];
% A=[1,1;-1,-1;-1,0];
% b=[2;1;1/2];
% [q,r]=qr(A);
% x0=[-2/3;2/3];
% [xk,r0,rk,fk,fm,fr]=FM(x0,q,r,A,b);
% [xk1,r0,rk,fk,fm,fr]=FM(xk,q,r,A,b);
% [xk2,r0,rk,fk,fm,fr]=FM(xk1,q,r,A,b);
% [xk3,r0,rk,fk,fm,fr]=FM(xk2,q,r,A,b);
% [xk4,r0,rk,fk,fm,fr]=FM(xk3,q,r,A,b);
% xa=[xk,xk1,xk2,xk3,xk4];
% for i=1:30
%     [xk4,r0,rk,fk,fm,fr]=FM(xa(:,4+i),q,r,A,b);
%     xa=[xa,xk4];
% end
% xa
