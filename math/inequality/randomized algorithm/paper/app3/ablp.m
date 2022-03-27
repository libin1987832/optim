function [x,arr]=ablp(A1,A2,l,u,amax,bmax,da,db,debug)
[m1,n]=size(A1);
[m2,n]=size(A2);
c=zeros(m1+n,1);
c((n+1):m1+n,1)=ones((n+1):m1+n,1);
alpha0=0;
beta0=0;
alpha=alpha0;
beta=beta0;
A=[A1,-diag(u);-A2,zeros(m1,m1);zeros(1,n) c'];
while 1
b=[zeros(m1,1);-l;m1*(1+alpha*beta)];
lb=zeros(m1+n,1);
ub=[Inf*ones(n,1);(1+beta)*ones(m1,0)];
[x,fval,flag] = linprog(c,A,b,[],[],lb,ub);
if flag==1
    break;
end
r=beta+db;
if r<=bmax
   beta=r;
end
if r>bmax
    beta=beta0;
    
end