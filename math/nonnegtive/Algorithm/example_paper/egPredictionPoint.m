addpath('../FM')
addpath('../hybrid')
addpath('../util')
A=[1,1;-1,-1;-1,0;6,3];
b=[1;1;0.5;-2];
[sm,sn]=size(A);
[Q,r]=qr(A);
Qn=Q(:,1:sn);
QQn=Qn*Qn';
nIter=3;
QQ=nIter*QQn;
A=[1,1;-1,-1;-1,0;-6,-3];
b=[1;1;0.5;2];
% output
d=lineData(A,b,[-2,1],[-1.5,3]);
line(d(:,1:2)',d(:,3:4)')
line([-2/3,1],[2/3,-1],'LineStyle','--');
hold on
xx1=[];
yy1=[];
xx2=[];
yy2=[];
xxx=[];
%  for x=-2:0.1:1
%      for y=-1.5:0.2:3
x=-0.6;
y=0.7;
x0=[x;y];
        fm=b-A*x0;
        ssign=getBn(QQ,fm,I);
        if ssign==sm
            xx1=[xx1;x];
            yy1=[yy1;y];
        else
            xx2=[xx2;x];
            yy2=[yy2;y];
        end
        [optimaln,real]=testReason(x0,A,b,nIter);
        if optimaln==nIter
            xxx=[xxx [x0;5]];
        elseif real==nIter
            xxx=[xxx [x0;optimaln]];
        else
            xxx=[xxx [x0;-1]];
        end
%     end
% end
% express neton's method good
plot(xx1,yy1,'r*');
hold on 
% express FM method good
plot(xx2,yy2,'bo');





