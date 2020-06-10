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
[xk2,err,indexG,indexN]=hybridorigin(x0,nmax,nf,C,q);
count=indexG/nf;
[xpa,errpa,indexpa,indexNpa]=PA(x0,nmax,nf,C,q);
countpa=indexpa/nf;
 xkA= splitS_our(C,q,1,x0,indexG);
 errsp=test_valid(C,q,xkA(:,indexG));
%  errspa=test_valid(C,q,all(:,indexG));
errspaa=test_valid(C,q,xk2);
% c2=cell2mat(err(4,2));
% s2=subspacesearch(c2,C,q);
% [e1,is1]=subspaceVaild(xs,C,q);
% [e2,is2]=subspaceVaild(s2,C,q);
% [e3,is3]=subspaceVaild(c2,C,q);
% norm(s2-xs)
% c1=xs;
% [n1,n2]=checkEq(c2,xs);
% c2(c2>0)=1;
% c1(c1>0)=1;
% sum(abs(c2-c1))
