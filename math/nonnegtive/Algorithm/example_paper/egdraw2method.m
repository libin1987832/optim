% compare simple iterator and newton method 
% A=[1,1;-1,-1;1,0;6,3];
% b=[1;1;-0.5;-2];
addpath('../FM');
addpath('../util');
addpath('../hybrid');
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

for x=-2.5:0.1:0
    for y=-1.5:0.2:3
        
        % form a point
        x0=[x;y];
        % in the expose face 
        %x0=[-2.5;1.7]
        % not in
        %x0=[-2;0];
        xF=x0;
        % two step compare
        for j=1:2
        [xk,r0,rk,fkFM,fm,fr]=FM(xF,q,r,A,b);
        xF=xk;
        end
        [xk,rk,fkssqr,f0,lambe]=ssqr(x0,A,b);
        % if the reduction is more ,then add the point
    if  fkFM>fkssqr
        % newtow good
        xx1=[xx1, x0(1)];
        yy1=[yy1, x0(2)];
    else
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

