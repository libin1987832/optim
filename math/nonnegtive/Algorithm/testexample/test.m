% A=diag([10^15,1]);
% % tf = issymmetric(A);
% 
% d = eig(A);
% tol=length(d)*eps(max(d))
% isposdef = all(d > tol)
% issemidef = all(d > -tol)
% 
% 
% n=5;
% A=rand(n+2,n);
% [q,r]=qr(A);
% Q1=q(:,1:3);
% N=diag([1,1,0,0,1,1,1])
% diag(ones(7,1))-Q1*Q1'*N

addpath('../FM')
addpath('../hybrid')
addpath('../util')
%  A=[1,1;-1,-1;-1,0;1,0;6,3];
%  b=[1;1;0.5;-3;2];

%  A=[1,1;-1,-1;-1,0;1,0];
%  b=[1;1;0.5;-3];

A=[1,1;-1,-1;-1,0;1,0;6,3];
b=[1;1;0.5;-1.5;-1.5];

 x0=[0;-0.5];

 
d=lineData(A,b,[-2,1],[-1.5,3]);
line(d(:,1:2)',d(:,3:4)')
line([-2/3,1],[2/3,-1],'LineStyle','--');
hold on 
[Q,R]=qr(A);
q=Q(:,1:2);
qq=q*q';
zs=[0,0,0,1,0]';
qzs=qq*zs;
qb=b-qq*b;
fmV=[];
for i=1:20
%     c=num2str(i);
%  plot(x0(1),x0(2),'*')
%  c=[' ',c];
% 
% text(x0(1),x0(2),c)
%  hold on
 fm0=A*x0-b;
 fmV=[fmV,fm0]
 fm0(fm0<0)=0;
 [xk,r0,rk,fk,fm,fr]=FM(x0,Q,R,A,b);
 x0=xk
%  qq(3,5)*fm0(5)+qq(3,4)*(fm0(4)-zs(4))
%  qq(5,5)*fm0(5)+qq(5,4)*(fm0(4)-zs(4))
%  qq(3,4)*(fm0(4)-zs(4))
%  qq(4,3)*(fm0(4)-zs(4))+qq(4,)
end
fmV
B=[3,4,5];
q3B=qq(3,B);
qBB=qq(B,B);
 x00=[0;-0.5];
 zs0=A*x00-b;
 zs0(zs0<0)=0;
 zsB=zs0-zs;
 qzs=qBB*zsB(B);
 qzs(qzs<0)=0;
z33=q3B*qzs
z34=qq(4,B)*qzs+zs(4)
z35=qq(5,B)*qzs
z=qBB*qzs+qBB*zs(B)