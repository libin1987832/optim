n=1000;
M=randn(n,n)*10^3;
M=M-diag(diag(M));
colM=sum(abs(M),2);
M=M+diag(colM);
xs=randn(n,1);
I1=(xs<0);
xs(I1)=0;
I2=(xs>0);
q=rand(n,1);
qt=M*xs;
q(xs>0)=-qt(xs>0);
% degency so +0.1
q(xs==0)=max(abs(qt))+0.1;
nmax=500;
etc=0.5;
ete=2;
rou=0.99;
trmax=1e12;
trr=1;
x0=ones(n,1)*100;
[x0,iter,nss]=asysub(x0,M,q,nmax,etc,ete,trr,trmax,rou);

