clc
clear
addpath('../FM')
addpath('../IFM')
addpath('../strategies')
addpath('../hybrid')
addpath('../linearsolve')
maxIter=300;
nmax=500;
Arecord=[];
     bnf=5;
 enf=30;
 nnf=enf-bnf+1;
    for m=200:100:200
   %m=1000;
    for ratio=0.6:0.2:0.8
      %  for batch=1:10
            nf=3;
         %for nf=bnf:enf
        n=floor(ratio*m);
         A=2*rand(m,n)-1;
         b=2*rand(m,1)-1;
        x0=zeros(n,1);
%          load('testA');
        xr=rand(n,1);
%         t=clock;
%         [q,r]=qr(A);
%         tq=etime(clock,t);
%        t=clock;
%         r\(q(1:n,:)'*xr);
%         tr=etime(clock,t);
%         t=clock;
%         for i=1:10
%         A*xr;
%         end
%         t1=etime(clock,t);
%         t=clock;
%         pinv(A);
%         t2=etime(clock,t);
%         fprintf('time A*x:%g,pinv:%g,qr:%g,r*q*x:%g\n',t1,t2,tq,tr);
        xs=-1;
        
%         [xkR,xkR2,countFR,countNWR,bNWR,tfR,vkR,rkArrR]=residualR(x0,A,b,maxIter);
%          dR=norm(xkR-xs);
%          rkR=b-A*xkR;
%          rkR(rkR<0)=0;
%          rkR=norm(rkR);
%          gR=norm(A'*rkR);
%          dR2=norm(xkR2-xs);
%          rkR2=b-A*xkR2;
%          rkR2(rkR2<0)=0;
%          gR2=norm(A'*rkR2);
% %         fprintf('resdual$ %d \\times %d $ & %g & %g & %4.2f & %d & %d &\n',m,n,rkR,gR,tfR,countFR,countNWR);
         
%          [xkA,rkA,countFA,countNA,bNWA,tfA,vkA,Arr]=als(x0,A,b,maxIter);
%          dA=norm(xkA-xs);
%          rkA=b-A*xkA;
%          rkA(rkA<0)=0;
%          gA=norm(A'*rkA);
%          fprintf('ALS$ %d \\times %d $ & %g & %g & %4.2f & %d & %d & %d\n',m,n,dA,gA,tfA,countFA,countNA,Arr(1,end));
%         Arecord=[Arecord;m n Arr(1,end) gA];
        
%          [xkpa,rkpa,countFMpa,countNWpa,beginNWpa,tfpa,vkpa]=pina(x0,A,b,maxIter);
%          dpa=norm(xkpa-xs);
%          rkpa=b-A*xkpa;
%          rkpa(rkpa<0)=0;
%          gpa=norm(A'*rkpa);
%          fprintf('pina$ %d \\times %d $ & %g & %g & %4.2f & %d & %d &\n',m,n,dpa,gpa,tfpa,countFMpa,countNWpa);
         
%          [xkhan,rkhan,countFMhan,countNWhan,beginNWhan,tfhan,vkhan,rkArrh]=han(x0,A,b,maxIter);
%          dhanD=norm(xkhan-xs);
%          rkhan=b-A*xkhan;
%          rkhan(rkhan<0)=0;
%          dhan=norm(rkhan);
%          ghan=norm(A'*rkhan);
 %        fprintf('han$ %d \\times %d $ & %g & %g & %4.2f & %d & %d &\n',m,n,dhan,ghan,tfhan,countFMhan,countNWhan);
                  
%          [xkLei,rkLei,countFMLei,countNWLei,beginNWLei,tfLei,vkLei]=Lei(x0,A,b,maxIter);
%          dLeiD=norm(xkLei-xs);
%          rkLei=b-A*xkLei;
%          rkLei(rkLei<0)=0;
%          dLei=norm(rkLei);
%          gLei=norm(A'*rkLei);
  %       fprintf('Lei$ %d \\times %d $ & %g & %g & %4.2f & %d & %d &\n',m,n,dLei,gLei,tfLei,countFMLei,countNWLei);
         
        [xkD,rkD,countFD,countND,bNWD,tfD,vkD,rkArrD]=Dax(x0,A,b,maxIter);
          xs=-1;
         dDD=norm(xkD-xs);
         rkD=b-A*xkD;
         rkD(rkD<0)=0;
         dD=norm(rkD);
         gD=norm(A'*rkD);
         fprintf('Dax$ %d \\times %d $ & %g & %g & %4.2f & %d & %d & %d\n',m,n,dD,gD,tfD,countFD,countND,bNWD);
        type=1;
        [xkG,rkG,countFG,countNG,bNWG,tfG,vkG,rkArrG]=gradientFM_i(x0,A,b,nf,1e-8,maxIter,xs,type);
          dGD=norm(xkG-xs);
          rkG=b-A*xkG;
          rkG(rkG<0)=0;
          dG=norm(rkG);
          gG=norm(A'*rkG);
          fprintf('grad$ %d \\times %d $ & %g & %g & %4.2f & %g & %g & %g &\n',m,n,dG,gG,tfG,countFG,countNG,bNWG);
 
     [xkC,rkC,countFMC,countNWC,beginNWC,tfC,vkC,rkArrC]=contraction_i(x0,A,b,nf,0.8,maxIter,xs,type);
         dCD=norm(xkC-xs);
        rkC=b-A*xkC;
          rkC(rkC<0)=0;
         dC=norm(rkC);
          gC=norm(A'*rkC);
          fprintf('con$ %d \\times %d $ & %g & %g & %4.2f & %g & %g & %g &\n',m,n,dC,gC,tfC,countFMC,countNWC,beginNWC);

       [xkP,rkP,countFP,countNP,bNWP,tfP,vkP,rkArrP]=predictFM_i(x0,A,b,nf,3,maxIter,xs,type);
        dPD=norm(xkP-xs);
        rkP=b-A*xkP;
          rkP(rkP<0)=0;
        dP=norm(rkP);
          gP=norm(A'*rkP);
         fprintf('pred$ %d \\times %d $ & %g & %g & %4.2f & %g & %g & %g &\n',m,n,dP,gP,tfP,countFP,countNP,bNWP); 
%          record=[m,n,gR,tfR,countFR,countNWR,nf;...
%                  m,n,gA,tfA,countFA,countNA,nf;...
%                  m,n,gD,tfD,countFD,countND,nf;...
%                  m,n,gG,tfG,countFG,countNG,nf;...
%                  m,n,gC,tfC,countFMC,countNWC,nf;...
%                  m,n,gP,tfP,countFP,countNP,nf...
%                 ];
%         Arecord=[Arecord;record];
         fprintf('$ %d \\times %d $ & %g & %g & %g & %g & %g & %g & %g & %g\n',m,n,gD,tfD,gC,tfC,gG,tfG,gP,tfP); 
   %     fprintf('%d\\time %d & %g & %g & %g & %g & %g & %g & %g & %g &\n',m,n,gD,tfD,gC,tfC,gG,tfG,gP,tfP);
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
%     end
end
    end

subplot(2,3,1)
rkArrR=rkArrR(rkArrR>0);
[m,n]=size(rkArrR);
plot(1:m,rkArrR)
subplot(2,3,2)
rkArrh=rkArrh(rkArrh>0);
[m,n]=size(rkArrh');
plot(1:m,rkArrh)
subplot(2,3,3)
rkArrD=rkArrD(rkArrD>0);
[m,n]=size(rkArrD);
plot(1:m,rkArrD)

subplot(2,3,4)
rkArrC=rkArrC(rkArrC>0);
[m,n]=size(rkArrC);
plot(1:m,rkArrC)
subplot(2,3,5)
rkArrG=rkArrG(rkArrG>0);
[m,n]=size(rkArrG);
plot(1:m,rkArrG)

subplot(2,3,6)
rkArrP=rkArrP(rkArrP>0);
[m,n]=size(rkArrP);
plot(1:m,rkArrP)
    
    
    % save(['ff' num2str(batch) '.mat'],'Arecord')
% end
% clear
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