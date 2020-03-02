% paper rand [1,1] example

m=300;
ratio=0.1;
n=ratio*m;
A=2*rand(m,n)-1;
b=2*rand(m,1)-1;
x0=zeros(n,1);

[xk1,fk1,xkArr1,countF1,countN1]=hybrid1(x0,A,b);
[xk2,fk2,xkArr2,countF2,countN2]=hybrid2(x0,A,b);

% [xk2,fk2,xkArr2,countF2,countN2]=hybrid4(x0,A,b);

