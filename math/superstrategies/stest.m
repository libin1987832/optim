clc
clear
addpath('../FM')
addpath('../IFM')
addpath('../strategies')
addpath('../hybrid')
addpath('../linearsolve')
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

Arecord=[];
bnf=10;
 enf=30;
 nnf=enf-bnf+1;
 for m=200:400:1000
    for ratio=0.1:0.1:1
        nf =30;
        for batch=1:10
        n=floor(ratio*m);
         A=2*rand(m,n)-1;
         b=2*rand(m,1)-1;
        x0=zeros(n,1);
        xr=rand(n,1);
        xs=-1;
        
        [xkR,xkR2,countFR,countNWR,bNWR,tfR,vkR]=residualR(x0,A,b,maxIter);
         dR=norm(xkR-xs);
         rkR=b-A*xkR;
         rkR(rkR<0)=0;
         gR=norm(A'*rkR);
         dR2=norm(xkR2-xs);
         rkR2=b-A*xkR2;
         rkR2(rkR2<0)=0;
         gR2=norm(A'*rkR2);
         fprintf('resdual$ %d \\times %d $ & %g & %g & %4.2f & %d & %d &\n',m,n,gR2,gR,tfR,countFR,countNWR);
         
         [xkA,rkA,countFA,countNA,bNWA,tfA,vkA,Arr]=als(x0,A,b,maxIter);
         dA=norm(xkA-xs);
         rkA=b-A*xkA;
         rkA(rkA<0)=0;
         gA=norm(A'*rkA);
         fprintf('ALS$ %d \\times %d $ & %g & %g & %4.2f & %d & %d & %d\n',m,n,dA,gA,tfA,countFA,countNA,Arr(1,end));
      
                 
         [xkhan,rkhan,countFMhan,countNWhan,beginNWhan,tfhan,vkhan]=han(x0,A,b,maxIter);
         dhan=norm(xkhan-xs);
         rkhan=b-A*xkhan;
         rkhan(rkhan<0)=0;
         ghan=norm(A'*rkhan);
         fprintf('han$ %d \\times %d $ & %g & %g & %4.2f & %d & %d &\n',m,n,dhan,ghan,tfhan,countFMhan,countNWhan);
                  
         [xkLei,rkLei,countFMLei,countNWLei,beginNWLei,tfLei,vkLei]=Lei(x0,A,b,maxIter);
         dLei=norm(xkLei-xs);
         rkLei=b-A*xkLei;
         rkLei(rkLei<0)=0;
         gLei=norm(A'*rkLei);
         fprintf('Lei$ %d \\times %d $ & %g & %g & %4.2f & %d & %d &\n',m,n,dLei,gLei,tfLei,countFMLei,countNWLei);
         
        [xkD,rkD,countFD,countND,bNWD,tfD,vkD]=Dax(x0,A,b,maxIter);
          xs=-1;
         dD=norm(xkD-xs);
         rkD=b-A*xkD;
         rkD(rkD<0)=0;
         gD=norm(A'*rkD);
         fprintf('Dax$ %d \\times %d $ & %g & %g & %4.2f & %d & %d & %d\n',m,n,dD,gD,tfD,countFD,countND,bNWD);
        
        [xkG,rkG,countFG,countNG,bNWG,tfG,vkG]=sgradientFM_i(x0,A,b,nf,1e-8,maxIter,xs,2);
          dG=norm(xkG-xs);
          rkG=b-A*xkG;
          rkG(rkG<0)=0;
          gG=norm(A'*rkG);
          fprintf('grad$ %d \\times %d $ & %g & %g & %4.2f & %g & %g & %g &\n',m,n,dG,gG,tfG,countFG,countNG,bNWG);
 
     [xkC,rkC,countFMC,countNWC,beginNWC,tfC,vkC]=scontraction_i(x0,A,b,nf,0.8,maxIter,xs,2);
         dC=norm(xkC-xs);
        rkC=b-A*xkC;
          rkC(rkC<0)=0;
         gC=norm(A'*rkC);
          fprintf('con$ %d \\times %d $ & %g & %g & %4.2f & %g & %g & %g &\n',m,n,dC,gC,tfC,countFMC,countNWC,beginNWC);

       [xkP,rkP,countFP,countNP,bNWP,tfP,vkP]=spredictFM_i(x0,A,b,nf,10,maxIter,xs,2);
        dP=norm(xkP-xs);
        rkP=b-A*xkP;
          rkP(rkP<0)=0;
        gP=norm(A'*rkP);
         fprintf('pred$ %d \\times %d $ & %g & %g & %4.2f & %g & %g & %g &\n',m,n,dP,gP,tfP,countFP,countNP,bNWP); 
         record=[m,n,gR,tfR,countFR,countNWR,nf;...
                 m,n,gA,tfA,countFA,countNA,nf;...
                 m,n,gD,tfD,countFD,countND,nf;...
                 m,n,gG,tfG,countFG,countNG,nf;...
                 m,n,gC,tfC,countFMC,countNWC,nf;...
                 m,n,gP,tfP,countFP,countNP,nf...
                ];
        Arecord=[Arecord record];
end
end
save(['ff' num2str(batch) '.mat'],'Arecord')
end
clear
% f1=load('ff1.mat');
% f2=load('ff2.mat');
% f3=load('ff3.mat');
% f4=load('ff4.mat');
% f5=load('ff5.mat');
% f6=load('ff6.mat');
% f7=load('ff7.mat');
% f8=load('ff8.mat');
% f9=load('ff9.mat');
% f10=load('ff10.mat');
% summ=[f1.Arecord(:,1:3),f2.Arecord(:,3),f3.Arecord(:,3),f4.Arecord(:,3),...
%     f5.Arecord(:,3),f6.Arecord(:,3),f7.Arecord(:,3),f8.Arecord(:,3)...
%     ,f9.Arecord(:,3),f10.Arecord(:,3)];
% summ2=[summ mean(summ(:,3:12),2) std(summ(:,3:12),0,2)];
% subplot(1,3,1)
% plot(0.1:0.1:1,summ2(1:10,13))
% subplot(1,3,2)
% plot(0.1:0.1:1,summ2(11:20,13))
% subplot(1,3,3)
% plot(0.1:0.1:1,summ2(21:30,13))

 %       [xs,fk,xkArr,countFM,countNW,Q]=hybrid1(x0,A,b,maxIter);
       % xkArr
%         [xkM,fkM,xkArrM,countFM,countNM]=hybridMPLSQR(x0,A,b,1,0.00001,maxIter);
%         [xkS,fkS,countFMS,countNWS]=hybridSplit(x0,A,b,maxIter,20,5,etc,ete,trr,trmax,rou);
%         [xk6,fk6,xkArr6,countF6,countN6]=hybrid6(x0,A,b,5,20,maxIter);

 %        [xkD,rkD,countFD,countND,bNWD,tfD,vkD]=Dax_GS(x0,A,b,maxIter);


% recordG=Arecord(4:6:end,:);
% recordC=Arecord(5:6:end,:);
% recordP=Arecord(6:6:end,:);
% subplot(3,2,1)
% plot(recordG(1:nnf,7),recordG(1:nnf,4),'--r*','linewidth',2,'markersize',5)
% subplot(3,2,2)
% plot(recordG(1+nnf:nnf+nnf,7),recordG(1+nnf:nnf+nnf,4),'--b*','linewidth',2,'markersize',5)
% subplot(3,2,3)
% plot(recordC(1:nnf,7),recordC(1:nnf,4),'--r*','linewidth',2,'markersize',5)
% subplot(3,2,4)
% plot(recordC(1+nnf:nnf+nnf,7),recordC(1+nnf:nnf+nnf,4),'--b*','linewidth',2,'markersize',5)
% subplot(3,2,5)
% plot(recordP(1:nnf,7),recordP(1:nnf,4),'--r*','linewidth',2,'markersize',5)
% subplot(3,2,6)
% plot(recordP(1+nnf:nnf+nnf,7),recordP(1+nnf:nnf+nnf,4),'--b*','linewidth',2,'markersize',5)


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

%         A=[1,1;-1,-1;-1,0;-6,-3];
%         b=[1;1;0.5;2];
%         x0=[-5;0];



 %save('nfrand.mat','Arecord')