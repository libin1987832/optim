str = 'bcsstk02.mtx bcsstk04.mtx';
strindex = [1 12; 14 25];
iter = [2,1];
for i = 1:size(strindex,1)
[A,rows,cols,entries,rep,field,symm] = mmread(str(1,strindex(i,1):strindex(i,2)));
n=cols;
B=speye(n);
x1=sparse(ones(n,1)./n);
t=clock;
[x,crit, iters, nitBB, error]=SPL(A, B, x1, iter(i),  1e-8, 1e-6, 1e-9, 0);
tf=etime(clock,t);
lambda = (x' * A * x) / (x' * x);
w = A * x - lambda * x;
disp(['spl:lambda=' num2str(lambda)  ',iter=' num2str(iters) ',dualfeasible='...
    num2str(min(w)) ',nitBB=' num2str(nitBB) ',nLsyst=' num2str(iters*nitBB) ',Crit=' num2str(crit) ',Cpu=' num2str(tf)])
end

