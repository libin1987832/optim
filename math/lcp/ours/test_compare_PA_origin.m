clc
clear

% n=500;
% A=randn(n);
% A=A'*A;
% B=0.1*eye(n);
% C=A+B;
% xs=randn(n,1);
% I1=(xs<0);
% xs(I1)=0;
% I2=(xs>0);
% cs2=sum(I2);
% q=rand(n,1);
% qt=C*xs;
% q(xs>0)=-qt(xs>0);
% q(xs==0)=max(abs(qt))+0.1;
% save('fpi1','C','xs','q','n','cs2')

 load('fpiq')
x0=ones(n,1);
nmax=100;
nf=10;
[xk2,err,indexG,indexN,all]=hybridorigin(x0,nmax,nf,C,q);
count=indexG/nf;
[xpa,errpa,indexpa,indexNpa]=PA(x0,nmax,nf,C,q);
countpa=indexpa/nf;
 xkA= splitS_our(C,q,1,x0,indexG);
 errsp=test_valid(C,q,xkA(:,indexG));
%  errspa=test_valid(C,q,all(:,indexG));
errspaa=test_valid(C,q,xk2);