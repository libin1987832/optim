clc
clear
str = 'bcsstk02.mtx bcsstk04.mtx';
strindex = [1 12; 14 25];
iterstr = [98 30;984 973];
for i = 1:size(strindex,1) 
[A,rows,cols,entries,rep,field,symm] = mmread(str(1,strindex(i,1):strindex(i,2)));
n=cols;
B=speye(n);
x1=sparse(ones(n,1)./n);
t=clock;
[x,crit, iters, nitBB, error]=SPL(B, A, x1, iterstr(i,2),  1e-8, 1e-6, 1e-9, 0, 1);
% figure;
% semilogy(1:iters,error(1:iters));
tfsql=etime(clock,t);
lambdasql = (x' * A * x) / (x' * B * x);
wsql = - A * x + lambdasql * B * x;
disp(['spl:lambda=' num2str(lambdasql)  ',iter=' num2str(iters) ',dualfeasible='...
    num2str(min(wsql)) ',nitBB=' num2str(nitBB) ',nLsyst=' num2str(iters*nitBB) ',Crit=' num2str(crit) ',Cpu=' num2str(tfsql)])
t=clock;
[x, iter, italg1] = FqpEicp(A, B, x1, iterstr(i,1), 1000, 200, 1, 1e-16, 1e-10, 1e-9);
tf=etime(clock,t);
lambda = (x' * A * x) / (x' * B * x);
w = -A * x + lambda * B* x;
disp(['fqp:lambda=' num2str(lambda)  ',it=' num2str(iter) ',ItAlg1=' num2str(italg1) ',dualfeasible='...
    num2str(min(w)) ',Cpu=' num2str(tf)])
end