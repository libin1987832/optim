%% read data
% [A,rows,cols,entries,rep,field,symm]=mmread('../util/well1033.mtx');
%addpath('./dataInequality/');
%addpath('./algorithmInequality/')
clear
clc
% format shortEng
addpath('./bramley/ineq')
m = 1000;
n = 100;
rangeMax = 2;
rangeMin = -2;
count = 1;
num_alg = 4;
maxIter = 10000;
record=zeros(num_alg*count,5);
tol = 1e-2;
% tols = 1e-5;
tolsa = (1e-1).^(1:2:12);%[1e-5:1e-1:1e-13];
[mtol,ntol] = size(tolsa);
arrspeed=zeros(num_alg-1,ntol);
revarr=zeros(ntol*count*num_alg,maxIter+2);
A = 2 * rand(m , n)-1;
b = 2 * rand(m , 1)-1;
x0 = zeros(n , 1);
 save(['A.mat'],'A','b');
% load('A.mat','A','b')
for k = 1:ntol
    tols= tolsa(k);
%     e=randn(1, n);
%     e(e<tol)=e(e<tol)+tol;
%     [~,l]=min(e);
%     e(l)=tol;
%     E = diag(e); % 只要最大除最小等于1000即可
%     E(1,1) = 10;
%     U = orth(randn(m, m));
%     V = orth(randn(n, n));
%     A = U*[E;zeros(m-n,n)]*V';
%        A = 2 * rand(m , n)-1;
       
    for j = 1:count
        %  A = 2 * rand(m , n)-1;
%         b = 2 * rand(m , 1)-1;
%         x0 = zeros(n , 1);
    
        nf = 5;
       str = ['F','C','P','R'];
        str2 = {'FM','CHA','RHA','PHA'};

        iterA=zeros(num_alg,maxIter+2);
        maxIterA = 0;
        fprintA=zeros(1,8);
        steplengthOrk = 3;
        fprintf('\\hline \n \\multirow{9}{*}{$ %d\\times %d $}',m,n);
        for i=1:num_alg
            if i ==1
                type = str(i);
                [xkD,flag,relres,iter,resvec,arvec,itersm,tfD]=otherAlg(A,b,x0,maxIter,type,tols);
                iter =iter+1;
%                 arvec
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
            if i ==1
                fprintf('& %s & %g & %g & (%d,%d)  & %g \\\\\n',type,dD,gD,iter,sumiter,tfD);
                iterA(i,1:iter)=resvec(1,iter);
            else
                fprintf('& %s & %g & %g & (%d,%d)  & %g \\\\\n',[type,'HA'],dD,gD,iter*nf,sumiter,tfD);
                iterA(i,1:iter*(nf+1))=resvec;
            end
            
        end
        
%         revarr((k-1)*count*num_alg+1:(k-1)*count*num_alg+num_alg,:) = iterA;
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
    t810=records(2:num_alg,:,5);
   % t810s=permute(t810,[3,2,1]);
    speed = bsxfun(@rdivide,t17,t810);
   %speed = t17./t810;
    speed_avg = squeeze(sum(speed,2)/count);
    arrspeed(:,k) = speed_avg;
end
[tolsa;arrspeed]
semilogx(tolsa,arrspeed(1,:),'b*-',tolsa,arrspeed(2,:),'ro-',tolsa,arrspeed(3,:),'g+-')
set(gca,'XDir','reverse')
% 标题标注 
title(['The speed up comparison of ours and ' str2{1}]) 
% 坐标轴标注 
xlabel('the accuracy') 
ylabel('the speed up') 
legend([str2{1} '/CHA'],[str2{1} '/RHA'],[str2{1} '/PHA']) 
%    end % for ration
%    end % for m