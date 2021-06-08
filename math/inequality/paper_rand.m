%% read data 
% [A,rows,cols,entries,rep,field,symm]=mmread('../util/well1033.mtx');
% 2000X800 is good
addpath('./dataInequality/');
addpath('./algorithmInequality/');
addpath('./class')
clear
clc
m = 1000;
n = 600;
rangeMax = 2;
rangeMin = -2;
count = 10;
record=zeros(5*count,5);

param = parameter();
for j = 1:count
     [A,b,x0]=randInequality(m,n,2,1);
%    A = 2 * rand(m , n)-1;
%    b = 2 * rand(m , 1)-1;
%    x0 = zeros(n , 1);
%    save('test.mat','A','b','x0');
%     load('test.mat','A','b','x0');
   str = ['D','U','C','R','P'];
% for the solution so here
% debug = 1;
% if debug
%     [xkh,rkh,countFMh,countNWh,beginNWh,tfh,vkh,rkArrh]=han(x0,A,b,maxIter);
%     [rk, rkh, dh, gh] = residual(A,b,x0);
%     fprintf("active:%d,tf:%g",vkh,tfh);
% end
    iterA=size(5,param.maxIter+2);
    maxIterA = 0;
    fprintA=zeros(1,8);
    disp(num2str(j));
    fprintf('\\hline \n \\multirow{4}{*}{$ %d\\times %d $}',m,n);
% there are four methods to run
    for i=1:5
        type = [str(i) 'HA'];
        param.type = type;
        [xkD,flag,relres,iter,resvec,arvec,itersm,tfD]=hybridNqr(A,b,x0,param);
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
        record((j-1)*5+i,:)=[dD,gD,iter*param.nf,sumiter,tfD*2];
% print for tex      
        % fprintf('%s $ %d \\times %d $ & %g & %g & %g & %g & %g & %g\n',type,m,n,dD,gD,tfD,iter*nf,sumiter,beginN(1));
        fprintf('& %s & %g & %g & (%d,%d)  & %g \\\\\n',type,dD,gD,iter*param.nf,sumiter,tfD*2);
       iterA(i,1:iter+1)=resvec;
        if maxIterA < iter+1
            maxIterA = iter+1;
        end
    end
end
records = reshape(record,5,count,5);
recordp = permute(records,[1,3,2]);
meanp = squeeze(sum(recordp,3)/count);
 fprintf('\\hline \n \\multirow{5}{*}{$ %d\\times %d $}',m,n);
for i =1:5
    if i==2
       fprintf('& %s & %g & %g & (%d,%d)  & %g \\\\\n','DHA($ \mu=n_f $)',meanp(i,1),meanp(i,2),round(meanp(i,3)),round(meanp(i,4)),meanp(i,5));
    else
        fprintf('& %s & %g & %g & (%d,%d)  & %g \\\\\n',[str(i) 'HA'],meanp(i,1),meanp(i,2),round(meanp(i,3)),round(meanp(i,4)),meanp(i,5));
    end
end      

% type=['r','c','k','g'];
% typet=['+','o','v','s'];
% beginp = 1;
% figure
% %maxIterA = 70;
% h=semilogy(beginp:maxIterA,iterA(1,beginp:maxIterA),'bx');
% h.LineStyle = '--';
% hold on
% for i=2:4
%     h=semilogy(beginp:maxIterA,iterA(i,beginp:maxIterA),[type(i) typet(i)]);
%     h.LineStyle = '--';
% end
% legend('DHA','CHA','RHA','PHA');
% xlabel('Iteration Number');
% ylabel('the norm of the gradient');

%    end % for ration
%    end % for m