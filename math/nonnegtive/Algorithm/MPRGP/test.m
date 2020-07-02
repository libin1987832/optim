%% test split
clc
clear

n=10;
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
save('fpi','C','xs','q','n');

% ressvv=test_valid(C,q,xs);
% [symmetric,posdef]=isPosdef(C);
% [symmetric2,posdef2]=isPosdef(TC-CC);
% load('fpiq')
x0=ones(n,1);

nmax=100;

max_iter = 50;
tol_rel  = 0.0;
tol_abs  = 0.0;

xst=test_bnd(C,q,xs);

[xkb,s,iter,Aopt]=qp_bnd(C,q);
[w,xka,retcode] = LCPSolve(C,q);
xkm=MPRGP(C,-q,x0,1,0.5,1e-10,0);
% xsa=norm(xs-xka);
% [res1,fx1]=test_valid(C,q,xka);
% xsb=norm(xs-xkb);

xspsor=norm(xs-xkpsor);
xs2=norm(xs-xk2);
xspa=norm(xs-xkpa);
xsor=norm(xs-xkor);
[ress,fxs]=test_valid(C,q,xs);
[resss,fxss]=test_valid(C,q,xks);
% [res,fx]=test_valid(C,q,xkb);
  
[res2,fx2]=test_valid(C,q,xkpsor);
[res3,fx3]=test_valid(C,q,xk2);
[res4,fx4]=test_valid(C,q,xkpa);
[resor,fxor]=test_valid(C,q,xkor);

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
%% test splitS pgs
% load('fpi')
% x0=ones(n,1);
% max_iter = 1;
% tol_rel  = 0.0;
% tol_abs  = 0.0;
% % [xkpgs err iter flag convergence msg] =  pgs(C, q, x0, max_iter, tol_rel, tol_abs, false);
% [xks,ress]=splitS(C,q,1.4,x0,10);
% [xkd,resd]=splitD(C,q,x0,20);
% [xkpsor err iter flag convergence msg] = psor(C, q, x0, 1.4, max_iter, tol_rel, tol_abs, false)
%% test split
% [xkpgs err iter flag convergence msg] =  pgs(C, q, x0, max_iter, tol_rel, tol_abs, false);
% [xk2,err]=splitForlcp(x0,nmax,C,q);
%     [res,fx]=test_valid(C,q,x0);
%     [res1,fx1]=test_valid(C,q,xk2);
%     [res2,fx2]=test_valid(C,q,xkpgs);
%     [res3,fx3]=test_valid(C,q,xs);