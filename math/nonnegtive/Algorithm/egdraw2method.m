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
for x=-2:0.1:1
    for y=-1.5:0.2:3
        x0=[x;y];
        [xk,r0,rk,fkFM,fm,fr]=FM(x0,q,r,A,b);
        [xk,rk,fkssqr,f0,lambe]=ssqr(x0,A,b);
    if  fkFM>fkssqr*1.0001
        plot(x,y,'ro');
    else
        plot(x,y,'b+');
    end
    hold on
    end
end

