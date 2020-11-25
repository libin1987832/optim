%% test split
clc
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
save('fpi','C','xs','q','n');

% ressvv=test_valid(C,q,xs);
% [symmetric,posdef]=isPosdef(C);
% [symmetric2,posdef2]=isPosdef(TC-CC);
% load('fpiq')
x0=ones(n,1);
e=max(eig(C));

nmax=100;

max_iter = 50;
tol_rel  = 0.0;
tol_abs  = 0.0;

xst=test_bnd(C,q,xs);
fxs=0.5*xs'*C*xs+q'*xs
[xkb,s,iter,Aopt]=qp_bnd(C,q);
[w,xka,retcode] = LCPSolve(C,q);
xkm=MPRGP(C,-q,x0,1,1/e,1e-10,0,nmax);
norm(xs-xkb)
norm(xs-xka)
norm(xs-xkm)