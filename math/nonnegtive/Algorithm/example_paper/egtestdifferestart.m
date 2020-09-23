addpath('../FM')
addpath('../hybrid')
addpath('../util')
clear;
A=[1,1;-1,-1;-1,0;-6,-3];
b=[1;1;0.5;2];
[Q,R]=qr(A);
xkA=[];
xkP=[];
rkA=[];
x0=[1,20]';
d=lineData(A,b,[-2,1],[-1.5,3]);
line(d(:,1:2)',d(:,3:4)')
%line([-2/3,1],[2/3,-1],'LineStyle','--');
line([-2,1],[2,-1],'LineStyle','--');
% while norm(x0)>2
% [xk,r0,rk,fk,fm,fr]=FM(x0,Q,R,A,b);
% x0=xk;
% end
hold on
for i=0:10
[xk,r0,rk,fk,fm,fr]=FM(x0,Q,R,A,b);
c=num2str(i);
plot(x0(1),x0(2),'*')
c=[' ',c];

text(x0(1),x0(2),c)
hold on
x0=xk;
rkn=b-A*x0;
p=find(rkn>0);
o=zeros(4,1);
o(p)=1;
rkA=[rkA norm(rk)];
xkP=[xkP o];
xkA=[xkA xk];
end
% 
 xkA
 rkA
% xkP