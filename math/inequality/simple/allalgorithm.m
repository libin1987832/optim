%% read data 
% [A,rows,cols,entries,rep,field,symm]=mmread('../util/well1033.mtx');
%addpath('./dataInequality/');
%addpath('./algorithmInequality/')
addpath('./bramley/ineq')
m = 1000;
n = 100;
rangeMax = 2;
rangeMin = -2;
count = 1;
num_alg = 9;
record=zeros(num_alg*count,5);
for j = 1:count
   A = 2 * rand(m , n)-1;
   b = 2 * rand(m , 1)-1;
   x0 = zeros(n , 1);
   maxIter = 300;
   nf = 5;
   str = ['F','I','P','H','B','D','C','R','P'];
    str2 = {'FM','IFM','PC','HAN','BW','DHA','CHA','RHA','PHA'};
    iterA=size(num_alg,maxIter+2);
    maxIterA = 0;
    fprintA=zeros(1,8);
    steplengthOrk = 3;
    fprintf('\\hline \n \\multirow{9}{*}{$ %d\\times %d $}',m,n);
    for i=1:num_alg
        if i <6
            type = str(i);
            [xkD,flag,relres,iter,resvec,arvec,itersm,tfD]=otherAlg(A,b,x0,maxIter,type); 
        else
             type = str(i);
            [xkD,flag,relres,iter,resvec,arvec,itersm,tfD]=hybridA(A,b,x0,steplengthOrk,maxIter,nf,[type,'HA']);
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
        if i <6
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
 fprintf('\\hline \n \\multirow{9}{*}{$ %d\\times %d $}',m,n);
for i =1:4
str = ['D','C','R','P'];
fprintf('& %s & %g & %g & (%d,%d)  & %g \\\\\n',[str(i),'HA'],meanp(i,1),meanp(i,2),round(meanp(i,3)),round(meanp(i,4)),meanp(i,5));
end      
%    end % for ration
%    end % for m