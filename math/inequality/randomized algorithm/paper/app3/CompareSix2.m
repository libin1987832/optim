clear
clc
m1=10000;
m2=6000;
n=2000;
al=0.2;
be=0.3;
A1=rand(m1,n);
A2=rand(m2,n);
x=rand(n,1);
d1=A1*x;
d2=A2*x;
dn=floor(m1*al);
u=d1;
u(1:dn,:)=u(1:dn,:)-be*rand(dn,1).*u(1:dn,:);

l=d2;
H=[-A1;A2;eye(n)];
b=[-u;l;zeros(n,1)];
x0=zeros(n,1);

maxIter=500000;

tol = 1e-5;
xst=zeros(1,5);
t=clock;
%[x1,alpha,beta,flag]=ablp(A1,A2,u,d2,0.5,0.5,0.1,0.1);
x1=zeros(n,1);
xst(1)=etime(clock,t);
t=clock;
[x2,iter,error_k,iter_k,index_k] = GuassSeidelNE(H, b, x0,0.5 ,maxIter,tol,[],0);
xst(2)=etime(clock,t);
t=clock;
[x3,iter,error_k,iter_k,index_k] = simpleGuassSeidelNE(H, b, x0,2 ,maxIter,tol,[],0);
xst(3)=etime(clock,t);
t=clock;
[x4,iter,error_k,iter_k,index_k] = randGuassSeidelNE(H, b, x0,2 ,maxIter,tol,[],0);
xst(4)=etime(clock,t);
maxIter=1;
t=clock;
[x5]=Fluence(A1,A2,u,l,dn,10,maxIter,'quadprog');
xst(5)=etime(clock,t);

xs=[x1 x2 x3 x4 x5];
% [alpha,beta]
str=['P','C','U','R','N'];
for i = 1:5
r = -u + A1 * xs(:,i);
rnum1=sum(r>0);
r = l - A2 * xs(:,i);
rnum2=sum(r>0);
r=b-H*xs(:,i);
r(r<0)=0;
r_GS = norm(r);
rum3=sum(x<0);
fprintf('& %sCD &%d &%d  & %g & %g \\\\\n', str(i), rnum1,rnum2,r_GS, xst(i));
end

%  
% Untitled
% & PCD &0 &3000 &0 & 13750.2 & 0 \\
% & RCD &2000 &2119 &0 & 1429.75 & 62.647 \\
% & PCD &1547 &2467 &0 & 1098.91 & 14.575 \\
% & CCD &1547 &2467 &0 & 1098.91 & 13.412 \\
% & NCD &1319 &2694 &0 & 1190.08 & 67.275 \\

% & PCD &0 &4000  & 23266.1 & 0 \\
% & CCD &2795 &2649  & 2645.23 & 36.381 \\
% & UCD &2277 &3204  & 1777.45 & 26.927 \\
% & RCD &2276 &3204  & 1777.45 & 29.19 \\
% & NCD &1816 &3628  & 1982.11 & 87.937 \\
