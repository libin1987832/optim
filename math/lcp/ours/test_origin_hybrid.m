%% test split
clear
n=100;
A=randn(n);
A=A'*A;
B=0.1*eye(n);
C=A+B;
xs=randn(n,1);
I1=(xs<0);
xs(I1)=0;
I2=(xs>0);
cs2=sum(I2);
q=rand(n,1);
qt=C*xs;
q(xs>0)=-qt(xs>0);
q(xs==0)=max(abs(qt))+0.1;
save('fpi1','C','xs','q','n','cs2')
% load('fpi')
x0=ones(n,1);
nmax=10;
[xk2,err,indexG,indexN,all]=hybridorigin(x0,nmax,C,q);
% [xs';all']