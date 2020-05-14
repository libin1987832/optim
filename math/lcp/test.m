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

max_iter = 10;
tol_rel  = 0.0;
tol_abs  = 0.0;

[xkb,s,iter,Aopt]=qp_bnd(C,q);
[w,xka,retcode] = LCPSolve(C,q);
[xkd,resd]=splitD(C,q,x0,1000);
[xks,ress]=splitS(C,q,1,x0,2000);
[xkpgs err iter flag convergence msg] =  pgs(C, q, x0, max_iter, tol_rel, tol_abs, false);
[xkpsor err iter flag convergence msg] = psor(C, q, x0, 1.4, max_iter, tol_rel, tol_abs, false)
[xk2,err]=splitForlcp(x0,nmax,C,q);
xsb=norm(xs-xkb)
xsa=norm(xs-xka)
xsd=norm(xs-xkd)
xss=norm(xs-xks)
xspgs=norm(xs-xkpgs)
xspsor=norm(xs-xkpsor)
xs2=norm(xs-xk2)
 resn=test_valid(C,q,xs)
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