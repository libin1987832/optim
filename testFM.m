A=[1 ;-1 ;0 ];
b=[2;3;1];
x0=[1];
[q,r]=qr(A);
[x1,f1]=FM(x0,q,r,A,b)
dFM(A,b,x1)


