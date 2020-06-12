%% test split
clc
clear
addpath('./other')
addpath('./symsub')
addpath('../other')
addpath('../symsub')

n=1000;
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
% load('fpiq')
CC=triu(C,1);
TC=tril(C);
ressvvf=func(C,q,xs);
ressvv=test_valid(C,q,xs);
[symmetric,posdef]=isPosdef(C);
[symmetric2,posdef2]=isPosdef(TC-CC);

x0=ones(n,1);
nmax=100;

max_iter = 50;
tol_rel  = 0.0;
tol_abs  = 0.0;


% [xkb,s,iter,Aopt]=qp_bnd(C,q);
%[w,xka,retcode] = LCPSolve(C,q);
% xsa=norm(xs-xka);
% [res1,fx1]=test_valid(C,q,xka);
[xkpsor err iter flag convergence msg] = psor(C, q, x0, 1, max_iter, tol_rel, tol_abs, false);
[xks,ress]=splitS(C,q,1.4,x0,10);
nf=10;
[xk2,err,index2]=splitForlcp(x0,nmax,nf,C,q);
[xkpa,errpa,indexpa1,indexpa2]=PA(x0,nmax,nf,C,q);
[xkor,error,indexor,indexNor]=hybridorigin(x0,nmax,nf,C,q);
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
disp(['RFNP:(',num2str(index2*(nf+1)),' ',num2str(index2),') err:',num2str(res3)])
disp(['PA:(',num2str(indexpa1),' ',num2str(indexpa2),') err:',num2str(res4)])
disp(['our:(',num2str(indexor),' ',num2str(indexNor),') err:',num2str(resor)])
disp(['GS:(',num2str(max_iter),' ',num2str(0),') err:',num2str(res2)])