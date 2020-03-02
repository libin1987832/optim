A=[1,1;-1,-1;1,0;6,3];
b=[1;1;-0.5;-2];

q1=[1/sqrt(39);-1/sqrt(39);1/sqrt(39);6/sqrt(39)];
q2=[19;-19;-20;-3]/sqrt(1131);
Q=[q1,q2];

x=[-3;-1];
r0=b-A*x;
r0A=r0;
r0A(r0A>0)=1;
r0A(r0A<0)=0;
NK=diag(r0A)
IQ=diag([1,1,1,1])-Q*Q'*NK;


[q,r]=qr(A);

[xk,r0,rk,fkFM,fm,fr]=FM(x,q,r,A,b);

r0IQ=b-A*x;
rkIQ=IQ*r0IQ;
rkIQS=rkIQ;
rkIQS(rkIQ>0)=1;
rkIQS(rkIQ<0)=0;
NK2=diag(rkIQS)

rkIQ2=b-A*xk;