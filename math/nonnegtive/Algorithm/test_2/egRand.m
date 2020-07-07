clc
clear
addpath('../FM')
addpath('../IFM')
addpath('../strategies')

tfA1=[];
dfA1=[];
tfA6=[];
dfA6=[];
dim={};
maxIter=500;
nmax=500;
etc=0.5;
ete=2;
rou=0.99;
trmax=1e12;
trr=1;

for m=100:10:120
for ratio=0.1:0.1:0.3
        n=ceil(ratio*m);
        A=2*rand(m,n)-1;
        b=2*rand(m,1)-1;
        x0=zeros(n,1);
        [xs,fk,xkArr,countFM,countNW,Q]=hybrid1(x0,A,b,maxIter);
        [xkM,fkM,xkArrM,countFM,countNM]=hybridMPLSQR(x0,A,b,1,0.00001,maxIter);
       [xkS,fkS,countFMS,countNWS]=hybridSplit(x0,A,b,maxIter,20,5,etc,ete,trr,trmax,rou);
        [xk6,fk6,xkArr6,countF6,countN6]=hybrid6(x0,A,b,5,20,maxIter);
        [xkD,rkD,countFD,countND,bNWD,tfD,vkD]=Dax(x0,A,b,maxIter);
        [xkG,rkG,countFG,countNG,bNWG,tfG,vkG]=gradientFM(x0,A,b,1,0.00001,maxIter);
        [xkC,rkC,countFC,countNC,bNWC,tfC,vkC]=contraction(x0,A,b,maxIter,20,5,etc,ete,trr,trmax,rou);
        [xkP,rkP,countFP,countNP,bNWP,tfP,vkP]=predictFM(x0,A,b,5,20,maxIter);
        dD=norm(xkD-xs);
        dG=norm(xkG-xs);
        dC=norm(xkC-xs);
        dP=norm(xkP-xs);
        gD=norm(A'*rkD);
        gG=norm(A'*rkG);
        gC=norm(A'*rkC);
        gP=norm(A'*rkP);
        fprintf('Dax$ %d \\times %d $ & %g & %g & %4.2f &\n',m,n,dD,gD,tfD);
        fprintf('grad$ %d \\times %d $ & %g & %g & %4.2f &\n',m,n,dG,gG,tfG);
        fprintf('pred$ %d \\times %d $ & %g & %g & %4.2f &\n',m,n,dP,gP,tfP);
        fprintf('con$ %d \\times %d $ & %g & %g & %4.2f &\n',m,n,dC,gC,tfC);
        
    
        
%         r=b-A*xk1;
%         r(r<0)=0;
%         df1=norm(A'*r);
%         %fprintf('$ %d \\times %d $ & %4.2f & %g & %d & %d & %4.2f &',m,n,0.5*(r'*r),df1,countF1,countN1,tf1);
%         r=b-A*xk6;
%         r(r<0)=0;
%         df6=norm(A'*r);
%         fprintf('%4.2f & %g & %d & %d & %4.2f\\\\%% %4.2f %4.2f\n',0.5*(r'*r),df6,countF6,countN6,tf6,tf1/tf6,log10(df1/df6));
%         tfA1=[tfA1 tf1];
%         tfA6=[tfA6 tf6];
%         dfA1=[dfA1 df1];
%         dfA6=[dfA6 df6];
%         dim=[dim num2str(n)];
     end
end
% xais=1:size(tfA1,2);
% itxais=5;
% subplot(1,2,1)
% %plot(1:size(dfA1,2),dfA1,'--o',1:size(dfA1,2),dfA6,'-o');
% semilogy(1:size(dfA1,2),dfA1,'--o',1:size(dfA1,2),dfA6,'-o');
% set(gca,'XTick',xais(1:itxais:end));
% set(gca,'XTicklabel',dim(1:itxais:end));% X��ļǺ�
% legend("Dax","Ours")
% xlabel('varibles');
% ylabel('the log of the norm of gradient ');
% title('the accuracy as the varibles increase');
% subplot(1,2,2)
% plot(1:size(tfA1,2),tfA1,'--o',1:size(tfA6,2),tfA6,'-o');
% set(gca,'XTick',xais(1:itxais:end));
% set(gca,'XTicklabel',dim(1:itxais:end));% X��ļǺ�
% legend("Dax","Ours")
% xlabel('varibles');
% ylabel('the time(s)');
% title('the cost time as the varibles increase');
% save('data6','dfA1','dfA6','tfA1','tfA6')