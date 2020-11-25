addpath('FM');
addpath('util');
addpath('hybrid');
A=[1,1;-1,-1;-1,0;-6,-3];
b=[1;1;0.5;2];

[q,r]=qr(A);
% output
d=lineData(A,b,[-2.5,0],[-1.5,3]);
line(d(:,1:2)',d(:,3:4)')

line([-2/3,-2.5],[2/3,2.5],'LineStyle','--');
hold on
xx1=[];
yy1=[];
xx2=[];
yy2=[];
nIter=3;
count=0;
for x=-2.5:0.1:0
    for y=-1.5:0.2:3
       count=count+11;
        x0=[x;y];  
         % in the expose face 
  %      x0=[-2.5;1.7]
        % not in
        %x0=[-2;0];
   %     x0=[-2.5;0.5]
    [optimaln,real]=testReason(x0,A,b,nIter);

        % if the reduction is more ,then add the point
    if  nIter==optimaln
        xx1=[xx1, x0(1)];
        yy1=[yy1, x0(2)];
    elseif nIter==real;
        xx2=[xx2,x0(1)];
        yy2=[yy2,x0(2)];
    end
    end
end
% express neton's method good
plot(xx1,yy1,'r*');
hold on 
% express FM method good
plot(xx2,yy2,'bo');


