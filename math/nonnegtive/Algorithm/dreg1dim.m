
A=[1;-1;1];
b=[1;1;3 ];
x0=5;
[Q,R]=qr(A);
x=-2:0.1:10 ;
r=b*ones(size(x))-A*x;
r(r<0)=0;
fk=0.5*(sum(r.*r,1));
plot(x,fk);
hold on
for i=1:10
    [xk,fk]=FM1(x0,Q,R,A,b);
    plot(xk,fk,'*');
    f0=fk;
    x0=xk;
    z=A*xk-b;
    z(z<0)=0;
    bz=b+z;  
    r=A*x-bz*ones(size(x));
    fk=0.5*(sum(r.*r,1));
    plot(x,fk);
    hold on
end
