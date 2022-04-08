m=60;
n=5;
% A = 2*rand(m,n)-1;
% b = 2*rand(m,1)-1;
 load('test1')
x00 = zeros(n,1);
maxiter = 300;
normr=zeros(6,maxiter);
normrow=[];
index=[];
 for i = 1:n
   normrow = [normrow,A(:,i)'*A(:,i)];
   index = [index,i];
end
weight = normrow/sum(normrow);
x1=x00;
Arr=[];
for k = 1:maxiter
    i = randsample(index,1,true,weight);
    Arr=[Arr,i];
    r=b-A*x1;
    r(r<0)=0;
    normr(1,k)=r'*r;
    x1(i)=x1(i)+lambda*A(:,i)'*r/(A(:,i)'*A(:,i));
end
Ar1=norm(A'*r)
N1=diag(r>0);
sr1=sum(r>0)
Ie=eye(m);
skpp=computerMSE(A,N1,lambda,0);
x1=x00;
for k = 1:maxiter
    i=Arr(k);
    r=b-A*x1;
    r(r<0)=0;
    if normr(3,k)==1
          rn(rn<0)=0;
        normr(2,k) = rn;
    else
        normr(2,k)=r'*r;
    end
    Nt=diag(r>0);
    if sum(Nt~=N1) == 0 & k+1<maxiter
        rr=r*r';
        rn=Ie(:)'*skpp*rr(:);
        normr(3,k+1)=1;
    end
    x1(i)=x1(i)+lambda*A(:,i)'*r/(A(:,i)'*A(:,i));
end

x2=x00;
for k = 1:maxiter
    i = mod(k,n) + 1; 
    r=b-A*x2;
    r(r<0)=0;
    normr(4,k)=r'*r;
    x2(i)=x2(i)+lambda*A(:,i)'*r/(A(:,i)'*A(:,i));
end
Ar2=norm(A'*r)
sr2=sum(r>0)
N2=diag(r>0);
sum(diag(N2)~=diag(N1))
skpc=computerCycle(A,N2,lambda,1);
x2=x00;
for k = 1:maxiter
   i = mod(k,n) + 1; 
    r=b-A*x2;
    r(r<0)=0;
    if normr(6,k)==1 & i == 1 & k+n<maxiter
        rn(rn<0)=0;
        normr(5,k) = rn'*rn;
    else
        normr(5,k)=r'*r;
    end
    Nt=diag(r>0);
    if sum(Nt~=N2) == 0 & i == 1 & k+n<maxiter
        r=b-A*x2;
        rn=skpc*r;
        normr(6,k+n)=1;
    end
    x2(i)=x2(i)+lambda*A(:,i)'*r/(A(:,i)'*A(:,i));
end
semilogy(1:maxiter,normr(1,:),'r',1:maxiter,normr(4,:),'b')