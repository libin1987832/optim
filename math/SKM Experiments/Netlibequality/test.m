clc
clear
load('lp_agg.mat'); 
n = size(C,2);
m = size(C,1);
C = Problem.A;
d = Problem.b;
error_x = inv(C'*C) * C' * d;
ep = randn(n,1);
C = [full(C);eye(n)];
d = [full(d);error_x+ep];
xinit = zeros(n,1);

error_x = inv(C'*C) * C' * d;
error  = C*error_x-d;
maxIter = 1e2;
tstart=cputime;
skmeq(C,d,xinit,4*norm(error,Inf),maxIter);
cputime-tstart
tstart=cputime;
randkaczmarz(C,d,xinit,4*norm(error,Inf),maxIter);
cputime-tstart


