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
count = 1;
num_alg = 4;
maxIter = 10000;
tol = 1e-2;
tols = 1e-12;
% tolsa = (1e-1).^(0:2:13);%[1e-5:1e-1:1e-13];
nfA = (3:1:30);
[mtol,ntol] = size(nfA);
record=zeros(num_alg*count*ntol,5);
arrspeed=zeros(num_alg-1,ntol);
revarr=zeros(ntol*count*num_alg,maxIter+2);
 A = 2 * rand(m , n)-1;
        b = 2 * rand(m , 1)-1;
        x0 = zeros(n , 1);
   %           save(['Anf.mat'],'A','b');
 load(['Anf4.mat']);
%load(['A3000_300.mat']);
likehoodNewton=zeros(4,ntol);
%           [xkD,flag,relres,iter,resvec,arvec,itersm,tfH] =otherAlg(A,b,x0,maxIter,'H',tols);
%           tfH
for k=1:ntol
nf= nfA(k);
%     e=randn(1, n);
%     e(e<tol)=e(e<tol)+tol;
%     [~,l]=min(e);
%     e(l)=tol;
%     E = diag(e); % 只要最大除最小等于1000即可
%     E(1,1) = 10;
%     U = orth(randn(m, m));
%     V = orth(randn(n, n));
%     A = U*[E;zeros(m-n,n)]*V'; 
    for j = 1:count
        %  A = 2 * rand(m , n)-1;
 
   
        str = ['U','C','R','P','H','H','B','D','U','C','R','P'];
        str2 = {'DHM($\mu=n_f$)','CHA','RHA','PHA','PC','HAN','BW','DHA','UHA','CHA','RHA','PHA'};

        iterA=zeros(num_alg,maxIter+2);
        maxIterA = 0;
        fprintA=zeros(1,8);
        steplengthOrk = 3;
nf
        fprintf('\\hline \n \\multirow{9}{*}{$ %d\\times %d $}',m,n);
        for i=1:num_alg
            type = str(i);
            [xkD,flag,relres,iter,resvec,arvec,itersm,tfD]=hybridA(A,b,x0,maxIter,nf,[type,'HA'],tols);
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
            record((k-1)*count*num_alg+(j-1)*num_alg+i,:)=[dD,gD,iter*nf,sumiter,tfD];
            % print for tex
           fprintf('& %s & %g & %g & (%d,%d)  & %g \\\\\n',[type,'HA'],dD,gD,iter*nf,sumiter,tfD);          
        end
    end
    records = reshape(record,num_alg,ntol,5);
%     meanp = squeeze(sum(recordp,2)/count);
%     fprintf('\\hline \n \\multirow{10}{*}{$ %d\\times %d $}',m,n);
%     for i =1:num_alg
%         fprintf('& %s & %g & %g & (%d,%d)  & %g \\\\\n',str2{i},meanp(i,1),meanp(i,2),round(meanp(i,3)),round(meanp(i,4)),meanp(i,5));
%     end
    
%     t17=records(1,:,5);
%     %t173=repmat(t17,1,1,3);
%     t810=records(2:num_alg,:,5);
%    % t810s=permute(t810,[3,2,1]);
%     speed = bsxfun(@rdivide,t17,t810);
%    %speed = t17./t810;
%     speed_avg = squeeze(sum(speed,2)/count);
%     arrspeed(:,k) = speed_avg;
    
%     likehoodNewton(:,k) =squeeze((records(:,k,4).*nf)./records(:,k,3));
end
format long

figure
% plot([0 nfA],[tfH,records(1,:,4)],'k.-',[0 nfA],[tfH.records(2,:,4)],'b*-',[0 nfA],[tfH records(3,:,4)],'ro-',[0 nfA],[tfH records(4,:,4)],'g+-')
% % % set(gca,'XTickLabel',tolsa);
% % % 标题标注
% set(gca,'YLim',[0 6]);%X轴的数据显示范围
%  title('The overall number of Newton type iterations with increasing nf') 
% % % 坐标轴标注 
% xlabel('nf') 
% ylabel('the overall number of Newton type algorithm') 
% legend('DHA(μ= nf)','CHA','RHA','PHA') 

figure
tfH =1.1
  plot([0 2 nfA],[1.1 0.67 records(1,:,5)],'k.-',[0 nfA],[tfH records(2,:,5)],'b*-',[0 2 nfA],[tfH 0.35 records(3,:,5)],'ro-',[0 2 nfA],[tfH 0.328 records(4,:,5)],'g+-')
% plot(nfA,records(1,:,5),'k.-',nfA,records(2,:,5),'b*-',nfA,records(3,:,5),'ro-',nfA,records(4,:,5),'g+-')
% % set(gca,'XTickLabel',tolsa);
% % 标题标注
 set(gca,'YLim',[0.2  1.1]);%X轴的数据显示范围

 %title('The performance of hybrid algorithms with increasing n_f') 
% % 坐标轴标注 
xlabel('n_f') 
ylabel('CPU(s)') 
legend('DHA(μ= n_f)','CHA','RHA','PHA') 

% [nfA;likehoodNewton]
% plot(nfA,likehoodNewton(1,:),'k.-',nfA,likehoodNewton(2,:),'b*-',nfA,likehoodNewton(3,:),'ro-',nfA,likehoodNewton(4,:),'g+-')
% % % set(gca,'XTickLabel',tolsa);
% % % 标题标注
% set(gca,'YLim',[0 1.5]);%X轴的数据显示范围
%  title('The likehood Newton of hybrid algorithm with increasing nf') 
% % % 坐标轴标注 
% xlabel('nf') 
% ylabel('the likehood of Newton type algorithm') 
% legend('DHA(μ= nf)','CHA','RHA','PHA') 

%  [nfA;arrspeed]
% plot(nfA,arrspeed(1,:),'b*-',nfA,arrspeed(2,:),'ro-',nfA,arrspeed(3,:),'g+-')
% % set(gca,'XTickLabel',tolsa);
% % 标题标注 
% title('Compared with other algorithms, the speed up of CHA is relation with the accuracy') 
% % 坐标轴标注 
% xlabel('the accuracy') 
% ylabel('the speed up') 
% legend('CHA','RHA','PHA') 

%    end % for ration
%    end % for m