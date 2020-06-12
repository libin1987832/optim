addpath('../ours')
addpath('../ours/predict')
addpath('../util')
clc
clear
n=100;
m=50;
A=randn(n);
A=A'*A;
B=0.5*eye(n);
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
x0r=repmat(xs,1,m);
x0r(I2,:)=x0r(I2,:)+rand(cs2,m)+ones(cs2,m);
x1r=rand(n,m);
x1r(x1r<0)=0;
save('fpi','C','xs','q','n','cs2','x0r','m','x1r');
%   load('fpi')
B=tril(C);
Ct=triu(C,1);
lambda=max(abs(eig(inv(B)*Ct)));

corret=[0,0];
wrong=[0,0];
tol_rel  = 1e-5;
tol_abs  = 1e-10;

data=zeros(m,2*n+2);
data(1,:)=[xs>0;xs;1;0]';

for j=2:2
    for i=1:m
        if j==1
            x0=x0r(:,i);
        else
            x0=x1r(:,i);
        end
        y2=checkEqS(xs,x0);
        [ynf0 err iter flag convergence msg] =  pgs(C,q, x0, 10, tol_rel, tol_abs, true);
        y1=checkEqS(ynf0,x0);
        [ynf1 err iter flag convergence msg] =  pgs(C,q, ynf0, 1, tol_rel, tol_abs, true);
        y2=checkEqS(ynf1,x0);
        [ynf2 err iter flag convergence msg] =  pgs(C,q, ynf1, 1, tol_rel, tol_abs, true);
        y2=checkEqS(ynf2,x0);
        xkA=[ynf0 ynf1 ynf2];
        [t,tlambda]=predict3(xkA(:,1),xkA(:,2),xkA(:,3));
        y1=checkEqS(t,x0);
        y2=checkEqS(xs,x0);
        data(i+1+m*j-m,:)=[t' x0' y1==y2 norm(ynf0-xs)];
        if y1==y2
            corret(j)=corret(j)+1;
        else
            wrong(j)=wrong(j)+1;
        end
    end
end
corret
[tlambda,lambda]
data;
iter=50;
Arr=zeros(n+2,iter+1);
x0=x0r(:,1);
Arr(:,1)=[x0r(:,1);test_valid(C,q,x0);0.5*x0'*C*x0+q'*x0];
for k=1:iter
    [ynf0 err iter flag convergence msg] =  pgs(C,q, x0, 1, tol_rel, tol_abs, true);
    x0=ynf0;
Arr(:,k+1)=[ynf0;err;0.5*ynf0'*C*ynf0+q'*ynf0];
end
Arr;
[s,z]=isPosdef(B-C)