%% test split
clear
% n=100;
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
%save('fpi1','C','xs','q','n','cs2')
load('fpi1')
x0=ones(n,1);
nmax=20;
nf=10;
[xk2,err,indexG,indexN,all]=hybridorigin(x0,nmax,nf,C,q);
count=indexG/nf;
[xk25,err5,indexG5,indexN5,all5]=hybridorigin(x0,nmax,5,C,q);
count5=indexG5/5;
cmax=0;
 [xkA,~] = splitS(C,q,1,x0,indexG*2,cmax);
% [xs';all']