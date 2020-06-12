addpath('../ours')
addpath('../ours/predict')
addpath('../util')
clc
clear
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
corret=0;
wrong=0;
tol_rel  = 1e-5;
tol_abs  = 1e-10;
for i=1:10
    x0=xs;
    x0(xs>0)=xs(xs>0)+rand(cs2,1)+ones(cs2,1);
    y2=checkEqS(xs,x0);
    [ynf0 err iter flag convergence msg] =  pgs(C,q, x0, 1, tol_rel, tol_abs, true);
    y1=checkEqS(ynf0,x0);
    [ynf1 err iter flag convergence msg] =  pgs(C,q, ynf0, 1, tol_rel, tol_abs, true);
    y2=checkEqS(ynf1,x0);
    [ynf2 err iter flag convergence msg] =  pgs(C,q, ynf1, 1, tol_rel, tol_abs, true);
    y2=checkEqS(ynf2,x0);
    xkA=[ynf0 ynf1 ynf2];
    t=predict2(xkA(:,1),xkA(:,2),xkA(:,3));
    y1=checkEqS(t,x0);
    y2=checkEqS(xs,x0);
    if y1==y2
        corret=corret+1;
    else
        wrong=wrong+1;
    end
end