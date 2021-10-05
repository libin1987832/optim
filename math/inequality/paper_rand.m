%% read data 
% [A,rows,cols,entries,rep,field,symm]=mmread('../util/well1033.mtx');
addpath('./dataInequality/');
addpath('./algorithmInequality/')
m = 1000;
n = 100;
rangeMax = 2;
rangeMin = -2;
count = 1;
record=zeros(4*count,5);
for j = 1:count
%     [A,b,x0]=randInequality(m,n,rangeMax,rangeMin);
   A = 2 * rand(m , n)-1;
   b = 2 * rand(m , 1)-1;
   x0 = zeros(n , 1);
   maxIter = 300;
   nf = 5;
   str = ['D','U','C','R','P'];
% for the solution so here
debug = 0;
if debug
    [xkh,rkh,countFMh,countNWh,beginNWh,tfh,vkh,rkArrh]=han(x0,A,b,maxIter);
    [rk, rkh, dh, gh] = residual(A,b,xkh);
    fprintf('active:%d ,%g',vkh,dh);
end
    iterA=size(4,maxIter+2);
    maxIterA = 0;
    fprintA=zeros(1,8);
    steplengthOrk = 3;
    fprintf('\\hline \n \\multirow{5}{*}{$ %d\\times %d $}',m,n);
% there are four methods to run
    for i=1:5
        type = str(i);
        [xkD,flag,relres,iter,resvec,arvec,itersm,tfD]=hybridA(A,b,x0,steplengthOrk,maxIter,nf,[type,'HA']);
  %      [rk, rkD, dD, gD] = residual(A,b,xkD);
        resvec = arvec;
% check the solution
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
        record((j-1)*5+i,:)=[dD,gD,iter*nf,sumiter,tfD];
% print for tex      
        % fprintf('%s $ %d \\times %d $ & %g & %g & %g & %g & %g & %g\n',type,m,n,dD,gD,tfD,iter*nf,sumiter,beginN(1));
        fprintf('& %s & %g & %g & (%d,%d)  & %g \\\\\n',[type,'HA'],dD,gD,iter*nf,sumiter,tfD);
       iterA(i,1:iter+1)=resvec;
        if maxIterA < iter+1
            maxIterA = iter+1;
        end
    end
end
records = reshape(record,5,count,5);
recordp = permute(records,[1,3,2]);
meanp = squeeze(sum(recordp,3)/count);
 fprintf('\\hline \n \\multirow{4}{*}{$ %d\\times %d $}',m,n);
for i =1:5
    str = ['D','U','C','R','P'];
fprintf('& %s & %g & %g & (%d,%d)  & %g \\\\\n',[str(i),'HA'],meanp(i,1),meanp(i,2),round(meanp(i,3)),round(meanp(i,4)),meanp(i,5));
end      
%    end % for ration
%    end % for m