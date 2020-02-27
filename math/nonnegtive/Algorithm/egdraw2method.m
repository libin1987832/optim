% compare simple iterator and newton method 
A=[1,1;-1,-1;1,0;6,3];
b=[1;1;-0.5;-2];

[q,r]=qr(A);

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

for x=-2:0.1:1
    for y=-1.5:0.2:3
        x0=[x;y];
        [xk,r0,rk,fkFM,fm,fr]=FM(x0,q,r,A,b);
        [xk,rk,fkssqr,f0,lambe]=ssqr(x0,A,b);
    if  fkFM>fkssqr*1.0001
        xx1=[xx1,x];
        yy1=[yy1,y];
    else
       xx2=[xx2,x];
        yy2=[yy2,y];
    end
    end
end
plot(xx1,yy1,'ro');
hold on 
plot(xx2,yy2,'b+');

