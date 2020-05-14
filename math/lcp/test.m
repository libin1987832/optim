%% test split
clear
n=10;
A=randn(n);
A=A'*A;
B=0.1*eye(n);
C=A+B;
xs=randn(n,1);
xs(xs<0)=0;
q=rand(n,1);
qt=C*xs;
q(xs>0)=-qt(xs>0);
q(xs==0)=qt(xs==0)+max(qt)+0.1;
save('fpi','C','xs','q','n')
% load('fpi')
x0=ones(n,1);
nmax=100;

[xkb,s,iter,Aopt]=qp_bnd(C,q);
[w,xka,retcode] = LCPSolve(C,q);
[xkd,resd]=splitD(C,q,x0,1000);
[xks,ress]=splitS(C,q,1,x0,2000);

xsb=norm(xs-xkb)
xsa=norm(xs-xka)
xsd=norm(xs-xkd)
xss=norm(xs-xks)
[xk2,err]=splitForlcp(x0,nmax,C,q);
xs2=norm(xs-xk2)
%% test qp_bnd
% n=3;
% A=randn(n)
% A=A'*A;
% B=0.1*eye(n);
% C=A+B;
% xs=randn(n,1);
% xs(xs<0)=0;
% q=randn(n,1);
% qt=-C*xs;
% q(xs>0)=qt(xs>0);
% q(xs==0)=qt(xs==0)+max(qt)+0.1;
%  save('fpi','C','xs','q')
% load('fpi')
% [xk,s,iter,Aopt]=qp_bnd(C,q);
% test_bnd(C,q,xs)
%% 