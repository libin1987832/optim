str = 'bcsstk02.mtx bcsstk04.mtx bcsstk19.mtx';
strindex = [1 12; 14 25; 27 38];
iter = [2000,100,200];
xa = [];
for i = 1:size(strindex,1)
[A,rows,cols,entries,rep,field,symm] = mmread(str(1,strindex(3,1):strindex(3,2)));
n=cols;
B=speye(n);
%x1=sparse(ones(n,1)./(n+1));
r=zeros(n,1);
for j=1:n
 r(j)=min(A(:,j)*B(j,j)-A(j,j)*B(:,j));
 if r(j)>=0
  disp('该矩阵不需要')
 end
end 
warning off 
s=find(r==max(r));
m=zeros(n,1);
m(s)=1;
x1=m;
t=clock;
[x,crit, iters, nitBB, error] = SPL(A, B, x1, iter(i),  1e-12, 1e-6, 1e-9, 0, 0);
tf=etime(clock,t);
xa = [xa, x];
%x = x +sprand(x)*1e-6;
lambda = (x' * A * x) / (x' * x);
w = A * x - lambda * x;
disp(['spl:lambda=' num2str(lambda)  ',iter=' num2str(iters) ',dualfeasible='...
    num2str(min(w)) ',nitBB=' num2str(nitBB) ',nLsyst=' num2str(iters*nitBB) ',Crit=' num2str(crit) ',Cpu=' num2str(tf)])
end
% full(xa)