addpath('./bramley/ineq')
m = 1000;
n = 400;
rangeMax = 2;
rangeMin = -2;
count = 1;
num_alg = 1;
record=zeros(num_alg*count,5);
tol = 1e-3;
% tols = 1e-5;
tolsa = (1e-1).^(0:9);%[1e-5:1e-1:1e-13];
[mtol,ntol] = size(tolsa);
arrspeed=zeros(num_alg-1,ntol);
format shortE
acc=zeros(ntol,100);
for k = 1:ntol
    tols= tolsa(k);
    e=randn(1, n);
    e(e<tol)=tol;
    E = diag(e); % 只要最大除最小等于1000即可
    E(1,1) = 10;
    U = orth(randn(m, m));
    V = orth(randn(n, n));
    A = U*[E;zeros(m-n,n)]*V';
    rank(A);
    for j = 1:count
        %  A = 2 * rand(m , n)-1;
        b = 2 * rand(m , 1)-1;
        x0 = zeros(n , 1);
        maxIter = 3000;
        nf = 5;
        str = ['H','C','H','C','F','H','B','D','U','C','R','P'];
        str2 = {'HAN','CHM','HAN','CHA','PC','HAN','BW','DHA','UHA','CHA','RHA','PHA'};
        iterA=size(num_alg,maxIter+2);
        maxIterA = 0;
        fprintA=zeros(1,8);
        steplengthOrk = 3;
  
                type = str(1);
                [xkD,flag,relres,iter,resvec,arvec,itersm,tfD]=otherAlg(A,b,x0,maxIter,type,tols);
               acc(k,1:size(arvec,2)+1)=[arvec tfD];
    end

end
