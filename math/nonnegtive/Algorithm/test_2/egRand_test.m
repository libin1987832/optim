clc
clear
addpath('../FM')
addpath('../IFM')
addpath('../strategies')
addpath('../hybrid')
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
trmax=1e2;
trr=1;

for m=1000:1000:3000
    for ratio=0.1:0.1:0.3
        n=ceil(ratio*m);
        A=2*rand(m,n)-1;
        b=2*rand(m,1)-1;
        x0=zeros(n,1);
 xs=-1;
 %       [xs,fk,xkArr,countFM,countNW,Q]=hybrid1(x0,A,b,maxIter);
       % xkArr
%         [xkM,fkM,xkArrM,countFM,countNM]=hybridMPLSQR(x0,A,b,1,0.00001,maxIter);
%         [xkS,fkS,countFMS,countNWS]=hybridSplit(x0,A,b,maxIter,20,5,etc,ete,trr,trmax,rou);
%         [xk6,fk6,xkArr6,countF6,countN6]=hybrid6(x0,A,b,5,20,maxIter);

       %  [xkD,rkD,countFD,countND,bNWD,tfD,vkD]=Dax(x0,A,b,maxIter);
%          dD=norm(xkD-xs);
%          gD=norm(A'*rkD);
%          fprintf('Dax$ %d \\times %d $ & %g & %g & %4.2f &\n',m,n,dD,gD,tfD);

       % [xkG,rkG,countFG,countNG,bNWG,tfG,vkG]=gradientFM(x0,A,b,1,0.00001,maxIter);
         [xkG,rkG,countFG,countNG,bNWG,tfG,vkG]=gradientFM_i(x0,A,b,1,0.99,1e-5,maxIter,-1);
          dG=norm(xkG-xs);
          rkG=b-A*xkG;
          rkG(rkG<0)=0;
          gG=norm(A'*rkG);
          fprintf('grad$ %d \\times %d $ & %g & %g & %4.2f & %g & %g & %g &\n',m,n,dG,gG,tfG,countFG,countNG,bNWG);

         [xkC,rkC,countFMC,countNWC,beginNWC,tfC,vkC]=contraction_i(x0,A,b,2,0.8,maxIter,-1);
%         [xkC,rkC,countFC,countNC,bNWC,tfC,vkC]=contraction_d(x0,A,b,maxIter,3,2,etc,ete,trr,trmax,rou);
   %      [xkC,rkC,countFC,countNC,bNWC,tfC,vkC]=contraction_d(x0,A,b,maxIter,20,5,etc,ete,trr,trmax,rou);
         dC=norm(xkC-xs);
        rkC=b-A*xkC;
          rkC(rkC<0)=0;
         gC=norm(A'*rkC);
          fprintf('con$ %d \\times %d $ & %g & %g & %4.2f & %g & %g & %g &\n',m,n,dC,gC,tfC,countFMC,countNWC,beginNWC);

  %      [xkPP,rkPP,countFPP,countNPP,bNWPP,tfPP,vkPP]=predictFM_d(x0,A,b,5,10,maxIter,xs);
        [xkP,rkP,countFP,countNP,bNWP,tfP,vkP]=predictFM(x0,A,b,5,10,maxIter);
        dP=norm(xkP-xs);
        rkP=b-A*xkP;
          rkP(rkP<0)=0;
        gP=norm(A'*rkP);
         fprintf('pred$ %d \\times %d $ & %g & %g & %4.2f & %g & %g & %g &\n',m,n,dP,gP,tfP,countFP,countNP,bNWP);        
    
        
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
% set(gca,'XTicklabel',dim(1:itxais:end));% X轴的记号
% legend("Dax","Ours")
% xlabel('varibles');
% ylabel('the log of the norm of gradient ');
% title('the accuracy as the varibles increase');
% subplot(1,2,2)
% plot(1:size(tfA1,2),tfA1,'--o',1:size(tfA6,2),tfA6,'-o');
% set(gca,'XTick',xais(1:itxais:end));
% set(gca,'XTicklabel',dim(1:itxais:end));% X轴的记号
% legend("Dax","Ours")
% xlabel('varibles');
% ylabel('the time(s)');
% title('the cost time as the varibles increase');
% save('data6','dfA1','dfA6','tfA1','tfA6')