addpath('../FM')
addpath('../hybrid')
addpath('../util')
clear;
A=[1,1;-1,-1;-1,0;-6,-3];
b=[1;1;0.5;2];
[Q,R]=qr(A);
QQ=Q(:,1:2)*Q(:,1:2)';
x=-0.6;
y=0.5;
x0=[x;y];
xkA=[];
xkP=[];
rkn=b-A*x0;
p=find(rkn>0);
o=zeros(4,1);
o(p)=1;
xkP=[xkP o];
Nk=rkn;
Nk(Nk>0)=1;
Nk(Nk<0)=0;
I=diag([1,1,1,1]);
Bn=I-QQ*diag(Nk);
r10=Bn^3*rkn
r102=(I-3*QQ*diag(Nk))*rkn
r1022=(I-2*QQ*diag(Nk)+1*QQ*diag(Nk)*QQ*diag(Nk))*rkn;
for i=0:2
[xk,r0,rk,fk,fm,fr]=FM(x0,Q,R,A,b);
% c=num2str(i);
% plot(x0(1),x0(2),'*')
% c=[' ',c];
% 
% text(x0(1),x0(2),c)
% hold on
x0=xk;

% p=find(rkn>0);
% o=zeros(4,1);
% o(p)=1;
% xkP=[xkP o];
% xkA=[xkA xk];
end
rkn=b-A*x0

d=lineData(A,b,[-2,1],[-1.5,3]);
line(d(:,1:2)',d(:,3:4)')
line([-2/3,1],[2/3,-1],'LineStyle','--');
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
xkP=[xkP o];
xkA=[xkA xk];
end
% 
% xkA
% xkP