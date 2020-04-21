addpath('../FM')
addpath('../hybrid')
addpath('../util')
A=[1,1;-1,-1;1,0;6,3];
b=[1;1;-0.5;-2];
[Q,r]=qr(A);
Qn=Q(:,1:2);
nIter=10;
QQ=nIter*Qn*Qn';
x=[-1;0.05;0];
y=[1;0.05:0];
A=[1,1;-1,-1;1,0;6,3];
b=[1;1;-0.5;-2];
% output
d=lineData(A,b,[-2,1],[-1.5,3]);
line(d(:,1:2)',d(:,3:4)')
line([-2/3,1],[2/3,-1],'LineStyle','--');
hold on
xx1=[];
yy1=[];
xx2=[];
yy2=[];
I=diag([1,1,1,1]);
for x=-2:0.1:1
    for y=-1.5:0.2:3
        x0=[x;y];
        fm=b-A*x0;
        ssign=getBn(QQ,fm,I);
        if ssign==4
            xx1=[xx1;x];
            yy1=[yy1;y];
        else
            xx2=[xx2;x];
            yy2=[yy2;y];
        end
    end
end
% express neton's method good
plot(xx1,yy1,'r*');
hold on 
% express FM method good
plot(xx2,yy2,'bo');





