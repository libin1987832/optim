addpath('../')
n=20;
A=randn(n);
A=A'*A;
B=0.1*eye(n);
C=A+B;
xs=randn(n,1);
xs(xs<0)=0;
q=rand(n,1);
qt=C*xs;
q(xs>0)=-qt(xs>0);
q(xs==0)=max(abs(qt))+0.1;
x0=ones(n,1);

xA=[];
I0=find(x0>0);
begin=1;
xAI=[];
iter=20;
I=zeros(1,iter);
xAI=[I];
xAI(begin,1)=1;
for i=2:iter
    [xks,ress]=splitS(C,q,1,x0,1);
    Iks=find(xks>0);
    e=setdiff(I0,Iks);
    if isempty(e)
        xAI(begin,i)=1;
    else
        xAI=[xAI;zeros(1,iter)];
        begin=begin+1;
        xAI(begin,i)=1;
        I0=Iks;
    end
    x0=xks;
    xA=[xA x0];
end
xAI
xA