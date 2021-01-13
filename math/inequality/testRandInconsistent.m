%% read data 
% [A,rows,cols,entries,rep,field,symm]=mmread('../util/well1033.mtx');
addpath('./dataInequality/');
addpath('./algorithmInequality/')
rand_example = false;
matrix_example = false;
seperat_example = false;

% [A,rows,cols,entries,rep,field,symm]=mmread('illc1850.mtx');
% [A,rows,cols,entries,rep,field,symm]=mmread('illc1033.mtx');
% [A,rows,cols,entries,rep,field,symm]=mmread('well1850.mtx');
% [A,rows,cols,entries,rep,field,symm]=mmread('well1033.mtx');
% m=rows;
% n=cols;
% % A(20:20:end,:)=0;
% %b=rand(rows,1);
% b=ones(m,1);
% b(1:2:end)=-1;

% [A1,b1,A2,b2,At1, At2, fm1, fm2, AL1, bL1, AL2, bL2] = readBreast(gamm1);
% [A1,b1,A2,b2,At1, At2, fm1, fm2, AL1, bL1, AL2, bL2] = readHeart(gamm1); 
% A = A1;
% b = b1;
% [m,n] = size(A);

% may be for many example 
% for m = 1000:1000:2000
%     for ratio = 0.6:0.1:0.8
    % m = 1000;
    % ratio=0.7;
%    n = floor( ratio * m);
%    A = 2 * rand(m , n)-1;
%    b = 2 * rand(m , 1)-1;

%% run program
% initialize the parameter
%    x0 = ones(n , 1);
m = 1000;
n = 200;
rangeMax = 2;
rangeMin = -1;
[A,b,x0]=randInequality(m,n,rangeMax,rangeMin);
    maxIter = 900;
    nf = 5;
    str = ['D','C','R','P'];
    ATA = A'*A;
    steplength = 1/(max(eig(ATA))+0.0001);
% for the solution so here    
    [xkh,rkh,countFMh,countNWh,beginNWh,tfh,vkh,rkArrh]=han(x0,A,b,maxIter);
    [rk, rkh, dh, gh] = residual(A,b,x0);
    %      rkh=b-A*xkh;
%      rkh(rkh<0)=0;
%     dh=norm(rkh);
%      gh=norm(A'*rkh);
  %  fprintf('han:%g %4.4f\n',gh,tfh);
%    [xk,rk,countFM,countNW,beginNW,tf,vk,Arr,rkrr]=als(x0,A,b,maxIter);
  %  fprintf('estiate value:%d\n',Arr(1,end));
    iterA=size(4,maxIter+2);
    maxIterA = 0;
    fprintA=zeros(1,8);
    record=zeros(4,5);
    steplengthOrk = 3;
    fprintf('\\hline \n \\multirow{4}{*}{$ %d\\times %d $}',m,n);
% there are four methods to run
    for i=1:4
        type = str(i);
%         
        [xkD,flag,relres,iter,resvec,arvec,itersm,tfD]=hybridNqr(A,b,x0,steplengthOrk,maxIter,nf,[type,'HA']);
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
        record(i,:)=[dD,gD,iter*nf,sumiter,tfD];
% print for tex      
        % fprintf('%s $ %d \\times %d $ & %g & %g & %g & %g & %g & %g\n',type,m,n,dD,gD,tfD,iter*nf,sumiter,beginN(1));
       fprintf('& %s & %g & %g & (%d,%d)  & %g \\\\\n',[type,'HA'],dD,gD,iter*nf,sumiter,tfD);
       iterA(i,1:iter+1)=resvec;
        if maxIterA < iter+1
            maxIterA = iter+1;
        end
    end
     
%    end % for ration
%    end % for m


[xkLei,rkLei,countFMLei,countNWLei,beginNWLei,tfLei,vkLei]=Lei(x0,A,b,maxIter,3);
% rkLei=b-A*xkLei;
% rkLei(rkLei<0)=0;
% dLei=norm(rkLei);
% gLei=norm(A'*rkLei);
[rk, rkLei,dLei, gLei] = residual(A,b,xkD);
fprintf('Lei$ %d \\times %d $ & %g & %g & %4.5f & %d & %d & %g,%g&\n',m,n,dLei,gLei,tfLei,countFMLei,countNWLei,sum(rkLei>1e-13),sum(rkLei>=0));

type=['r','c','k','g'];
typet=['+','o','v','s'];
beginp = 1;
figure
%maxIterA = 70;
h=semilogy(beginp:maxIterA,iterA(1,beginp:maxIterA),'bx');
h.LineStyle = '--';
hold on
for i=2:4
    h=semilogy(beginp:maxIterA,iterA(i,beginp:maxIterA),[type(i) typet(i)]);
    h.LineStyle = '--';
end
legend('DHA','GHA','RHA','PHA');
xlabel('Iteration Number');
ylabel('the norm of the gradient');
    