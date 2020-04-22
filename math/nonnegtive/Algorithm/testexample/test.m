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
%  b=[1;1;0.5;-3;-2];
 A=[1,1;-1,-1;-1,0;1,0];
 b=[1;1;0.5;-3];

 x0=[0;0.5];

 
d=lineData(A,b,[-2,1],[-1.5,3]);
line(d(:,1:2)',d(:,3:4)')
line([-2/3,1],[2/3,-1],'LineStyle','--');
hold on 
[Q,R]=qr(A);
q=Q(:,1:2);
qq=q*q';
zs=[0,0,0,2.5]';
qzs=qq*zs;
for i=1:10
 plot(x0(1),x0(2),'*')
 hold on
 [xk,r0,rk,fk,fm,fr]=FM(x0,Q,R,A,b);
 x0=xk
 qzs(:,4)-
end
