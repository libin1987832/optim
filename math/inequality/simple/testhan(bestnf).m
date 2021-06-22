addpath('./bramley/ineq')
clear
clc
m = 1000;
n = 200;
rangeMax = 2;
rangeMin = -2;
count = 1;
num_alg = 1;
record=zeros(num_alg*count,5);
tol = 1e-3;
% tols = 1e-5;
tolsa = (1e-1).^(12);%[1e-5:1e-1:1e-13];
[mtol,ntol] = size(tolsa);
arrspeed=zeros(num_alg-1,ntol);
format shortE
acc=zeros(ntol,100);
acc1=zeros(ntol,100);
inf_num=3;
anor= zeros(inf_num*ntol,300);

for k = 1:ntol
    tols= tolsa(k);
%     e=randn(1, n);
%     e(e<tol)=tol;
%     E = diag(e); % 只要最大除最小等于1000即可
%     E(1,1) = 10;
%     U = orth(randn(m, m));
%     V = orth(randn(n, n));
%     A = U*[E;zeros(m-n,n)]*V';
%     rank(A);
%        A = 2 * rand(m , n)-1;
%    b = 2 * rand(m , 1)-1;
load('A.mat')
[m,n]=size(A);
  x0 = zeros(n , 1);

    for j = 1:count
        %  A = 2 * rand(m , n)-1;
%         b = 2 * rand(m , 1)-1;
%         x0 = zeros(n , 1);
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
        rpkt = b-A*xkD;
        rpkt(rpkt<0)=0;
        It= rpkt>0;
        [Q,R]=qr(A);
       % [xk0o,flag,relres,iter,resvec,arvecN,itersm,tfD]=otherAlg(A,b,x0,maxIter,'N',tols);
%         acc1(k,1:size(arvecN,2)+1)=[arvecN tfD];
        for i = 1:200
            rk=b-A*x0;
            [xk,rk]=FMQR(x0,Q,R,A,b,rk);
            [xs,rpk,len,flag]=sm(A,b,n,rk,xk);

            I=rk>0;
            rka=rk;
            rka(rka<0)=0;
%             anor(k,i*3-2:i*3) = [norm(rpkt-rka) norm(xk-xkD) sum(abs(I-It))];
            rpk(rpk<0)=0;
%            anor(k,i) = [sum(abs(I-It))];
             anor(2*k-1:2*k,i) = [sum(abs(I-It));norm(A'*rpk)];
            x0=xk;
            if norm(xs-xkD)<1e-15
                break;
            end
        end
    end

end
