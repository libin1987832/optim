clc
clear
str = 'bcsstk02.mtx bcsstk04.mtx bcsstk19.mtx';
strindex = [1 12; 14 25; 27 38];
iter = [2000,100,200];
xa = [];
for i = 1:size(strindex,1)
[A,rows,cols,entries,rep,field,symm] = mmread(str(1,strindex(3,1):strindex(3,2)));
n = cols;
B = speye(n);
x1=sparse(ones(n,1)./(n+1));
M = B;
sigma0 = max(eigs(A));
dsigma = 0.01;
eps = 1e-6;
bisectEps = 1e-12;
maxIT = 100;
debug = 0;
t=clock;
[x, iter, error]  = sqpEicp(A, B, M, x1, sigma0, dsigma, eps, bisectEps, maxIT, debug);
tf=etime(clock,t);
lambda = (x' * A * x) / (x' * x);
w = A * x - lambda * x;
disp(['spl:lambda=' num2str(lambda)  ',iter=' num2str(iter) ',dualfeasible='...
    num2str(min(w)) ',Cpu=' num2str(tf)])
end
% full(xa)