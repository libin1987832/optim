%% read data
% [A,rows,cols,entries,rep,field,symm]=mmread('../util/well1033.mtx');
%addpath('./dataInequality/');
%addpath('./algorithmInequality/')
clear
clc
% format shortEng
addpath('./bramley/ineq')
m = 3000;
n = 300;
rangeMax = 2;
rangeMin = -2;

num_alg = 4;
maxIter = 10000;
record=zeros(num_alg, 6);
tol = 1e-2;

revarr=zeros(num_alg, maxIter+2);

e=randn(1, n);
e(e<tol)=e(e<tol)+tol;
[~,l]=min(e);
e(l)=tol;
E = diag(e); % 只要最大除最小等于1000即可
E(1,1) = 10;
U = orth(randn(m, m));
V = orth(randn(n, n));
% A = U*[E;zeros(m-n,n)]*V';

 A = 2 * rand(m , n)-1;
  b = 2 * rand(m , 1)-1;
x0 = zeros(n , 1);
% load('A2000_100.mat');
save('A3000_300.mat','A','b');
nf = 5;

str = ['D','U','C','R','P'];
str2 = {'DHA','UHA','CHA','RHA','PHA'};

iterA=zeros(num_alg,maxIter+2);
maxIterA = 0;
fprintA=zeros(1,8);
steplengthOrk = 3;
for i = 1:5 
    type = str(i);
    [xkD,flag,relres,iter,resvec,arvec,itersm,tfD]=hybridA(A,b,x0,maxIter,nf,[type,'HA'],1e-12);
    resvec = arvec;
    rkD=b-A*xkD;
    rkD(rkD<0)=0;
    dD=norm(rkD);
    gD=norm(A'*rkD);
    beginN=find(itersm>0);
    sumiter=sum(itersm>0);
    fprintA(2*i-1:2*i)=[gD,tfD];
    if isempty(beginN)
        beginN=0;
    end
    record(i,:)=[dD,gD,iter,iter*nf,sumiter,tfD];
    fprintf('& %s & %g & %g & (%d,%d,%d)  & %g \\\\\n',[type,'HA'],dD,gD,iter,iter*nf,sumiter,tfD);
end


