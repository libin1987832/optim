%% test split
clear
n=500;
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
save('fpi','C','xs','q','n')
%  load('fpi')
x0=ones(n,1);
nmax=100;

[xk2,err,index2]=splitForlcp(x0,nmax,C,q);
xs2=norm(xs-xk2)
[res2,fx3]=test_valid(C,q,xk2);
