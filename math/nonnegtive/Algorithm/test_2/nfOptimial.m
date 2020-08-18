clc
clear
addpath('../FM')
addpath('../IFM')
addpath('../strategies')
addpath('../hybrid')
addpath('../linearsolve')
tfA1=[];  



dim={};
maxIter=500;
nmax=500;
etc=0.5;
ete=2;
rou=0.99;
trmax=1e2;
trr=1;



Nrecord=[];
     bnf=10;
 enf=30;
 nnf=enf-bnf+1;
   % m=200:400:1000
   m=1000;
    %for ratio=0.1:0.1:1
    for n=[100,500,1000]
     for nf=bnf:enf
         Arecord=[];
        for batch=1:10
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
       %  Arecord=[Arecord;m n Arr(1,end) gA];
         
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
        
        [xkG,rkG,countFG,countNG,bNWG,tfG,vkG]=gradientFM_i(x0,A,b,nf,1e-8,maxIter,xs,2);
          dG=norm(xkG-xs);
          rkG=b-A*xkG;
          rkG(rkG<0)=0;
           gG=norm(A'*rkG);
          fprintf('grad$ %d \\times %d $ & %g & %g & %4.2f & %g & %g & %g &\n',m,n,dG,gG,tfG,countFG,countNG,bNWG);
 
     [xkC,rkC,countFMC,countNWC,beginNWC,tfC,vkC]=contraction_i(x0,A,b,nf,0.8,maxIter,xs,2);
         dC=norm(xkC-xs);
        rkC=b-A*xkC;
          rkC(rkC<0)=0;
         gC=norm(A'*rkC);
          fprintf('con$ %d \\times %d $ & %g & %g & %4.2f & %g & %g & %g &\n',m,n,dC,gC,tfC,countFMC,countNWC,beginNWC);

       [xkP,rkP,countFP,countNP,bNWP,tfP,vkP]=predictFM_i(x0,A,b,nf,10,maxIter,xs,2);
        dP=norm(xkP-xs);
        rkP=b-A*xkP;
          rkP(rkP<0)=0;
        gP=norm(A'*rkP);
         fprintf('pred$ %d \\times %d $ & %g & %g & %4.2f & %g & %g & %g &\n',m,n,dP,gP,tfP,countFP,countNP,bNWP);
         recordE=[gR,tfR,countFR,countNWR;...
                 gA,tfA,countFA,countNA;...
                 gD,tfD,countFD,countND;...
                 gG,tfG,countFG,countNG;...
                 gC,tfC,countFMC,countNWC;...
                 gP,tfP,countFP,countNP...
                ];
        Arecord=[Arecord recordE]; 
        end
        AArecord=[repmat([m,n,nf],6,1) Arecord];
     end
     Nrecord=[Nrecord;AArecord];
    end
tfsumme=Nrecord(:,5:4:end);
tfsme=[tfsumme mean(tfsumme,2) std(tfsumme,0,2)];
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

%         A=[1,1;-1,-1;-1,0;-6,-3];
%         b=[1;1;0.5;2];
%         x0=[-5;0];



 %save('nfrand.mat','Arecord')