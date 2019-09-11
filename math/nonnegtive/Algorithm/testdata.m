m=2000;
n=0.1*m;
A=rand(m,n)*2-ones(m,n);
b=rand(m,1)*2-ones(m,1);
x0=zeros(n,1);

[xk1,fk1,xkArr1]=hybrid1(x0,A,b);
[xk2,fk2,xkArr2]=hybrid2(x0,A,b);
[xk3,fk3,xkArr3]=hybrid3(x0,A,b);