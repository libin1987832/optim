%% read data
% [A,rows,cols,entries,rep,field,symm]=mmread('../util/well1033.mtx');
%addpath('./dataInequality/');
%addpath('./algorithmInequality/')
clear
clc
addpath('./bramley/ineq')
m = 1000;
n = 300;
rangeMax = 2;
rangeMin = -2;
count = 3;
num_alg = 2;
record=zeros(num_alg*count,5);
tol = 1e-5;
% tols = 1e-5;
tolsa = linspace(1e-6,1e-15,3);%[1e-5:1e-1:1e-13];
[mtol,ntol] = size(tolsa);
arrspeed=zeros(1,ntol);
for k = 1:ntol
    tols= tolsa(k);
    e=randn(1, n);
    e(e<tol)=tol;
    E = diag(e); % 只要最大除最小等于1000即可
    E(1,1) = 10;
    U = orth(randn(m, m));
    V = orth(randn(n, n));
    A = U*[E;zeros(m-n,n)]*V';
    rank(A)
    for j = 1:count
        %  A = 2 * rand(m , n)-1;
        b = 2 * rand(m , 1)-1;
        x0 = zeros(n , 1);
        maxIter = 500;
        nf = 5;
        str = ['F','C','F','H','B','D','U','C','R','P'];
        str2 = {'FM','CHA','PC','HAN','BW','DHA','UHA','CHA','RHA','PHA'};
        iterA=size(num_alg,maxIter+2);
        maxIterA = 0;
        fprintA=zeros(1,8);
        steplengthOrk = 3;
        fprintf('\\hline \n \\multirow{9}{*}{$ %d\\times %d $}',m,n);
        for i=1:num_alg
            if i == 1
                type = str(i);
                [xkD,flag,relres,iter,resvec,arvec,itersm,tfD]=otherAlg(A,b,x0,maxIter,type,tols);
            else
                type = str(i);
                [xkD,flag,relres,iter,resvec,arvec,itersm,tfD]=hybridA(A,b,x0,steplengthOrk,maxIter,nf,[type,'HA'],tols);
            end
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
            record((j-1)*num_alg+i,:)=[dD,gD,iter*nf,sumiter,tfD];
            % print for tex
            % fprintf('%s $ %d \\times %d $ & %g & %g & %g & %g & %g & %g\n',type,m,n,dD,gD,tfD,iter*nf,sumiter,beginN(1));
            if i == 1
                fprintf('& %s & %g & %g & (%d,%d)  & %g \\\\\n',type,dD,gD,iter*nf,sumiter,tfD);
            else
                fprintf('& %s & %g & %g & (%d,%d)  & %g \\\\\n',[type,'HA'],dD,gD,iter*nf,sumiter,tfD);
            end
            iterA(i,1:iter+1)=resvec;
        end
        
        
    end
    records = reshape(record,num_alg,count,5);
    recordp = permute(records,[1,3,2]);
    meanp = squeeze(sum(recordp,3)/count);
    fprintf('\\hline \n \\multirow{10}{*}{$ %d\\times %d $}',m,n);
    for i =1:num_alg
        fprintf('& %s & %g & %g & (%d,%d)  & %g \\\\\n',str2{i},meanp(i,1),meanp(i,2),round(meanp(i,3)),round(meanp(i,4)),meanp(i,5));
    end
    
    t17=records(1,:,5);
    %t173=repmat(t17,1,1,3);
    t810=records(2,:,5);
   % t810s=permute(t810,[3,2,1]);
   % speed = bsxfun(@rdivide,t17,t810s);
   speed = t17./t810;
    speed_avg = mean(speed);
    arrspeed(k) = speed_avg;
end
arrspeed
%    end % for ration
%    end % for m