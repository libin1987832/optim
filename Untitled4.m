A=[-1;-1;1];
b=[3;-3;2];
[q,r]=qr(A);
x0=-4
[x1,f1]=FM(x0,q,r,A,b)
dFM(A,b,x1)
% plotFM(A,b,x0,x1)