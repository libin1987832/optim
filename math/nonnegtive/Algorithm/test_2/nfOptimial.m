clc
clear
addpath('../FM')
addpath('../IFM')
addpath('../strategies')
addpath('../hybrid')
addpath('../linearsolve')
tfA1=[];  



dim={};
maxIter=200;



Nrecord=[];
     bnf=10;
 enf=30;
 nnf=enf-bnf+1;
   % m=200:400:1000
   m=1000;
    %for ratio=0.1:0.1:1
    for n=[100,500,1000]
     for nf=bnf: enf
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
        Nrecord=[Nrecord;AArecord];
     end
     
    end
tfsumme=Nrecord(:,5:4:end);
tfsme=[tfsumme mean(tfsumme,2) std(tfsumme,0,2)];

%%
Daxtf=tfsme(3:6:end,:);
gradienttf=tfsme(4:6:end,:);
constracttf=tfsme(5:6:end,:);
predicttf=tfsme(6:6:end,:);
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
 
 gradient:
 1000:[0.0457999999999996;0.0474000000000003;0.0471999999999995;0.0462999999999989;0.0475000000000005;0.0475000000000001;0.0473000000000011;0.0477000000000002;0.0476000000000011;0.0484000000000009;0.0494999999999996;0.0482000000000001;0.0471999999999997;0.0464999999999998;0.0473999999999997;0.0498999999999991;0.0466000000000003;0.0483000000000005;0.0481999999999997;0.0479000000000003;0.0484000000000000]
 500:[6.22320000000000;6.20790000000000;6.19490000000000;6.21320000000000;6.01930000000000;6.06250000000000;6.07490000000000;5.92710000000000;5.99220000000000;5.84390000000000;6.03860000000000;5.91530000000000;5.95950000000000;5.94740000000000;5.80260000000000;5.90090000000000;6.03760000000000;5.91870000000000;5.89810000000000;5.94410000000000;5.88250000000000]
 100:[0.0937000000000001;0.0951999999999991;0.0980999999999995;0.0936000000000007;0.104599999999999;0.113600000000000;0.0923000000000002;0.0881999999999991;0.0791000000000004;0.0896000000000008;0.112400000000000;0.112500000000001;0.102800000000000;0.104900000000000;0.120499999999998;0.0989000000000000;0.103700000000000;0.103200000000001;0.124300000000001;0.126600000000001;0.119000000000000]
 contract:
 1000[0.0457000000000006;0.0479999999999999;0.0470999999999989;0.0487999999999993;0.0484999999999997;0.0470999999999994;0.0483999999999996;0.0494000000000002;0.0458000000000000;0.0487999999999993;0.0492000000000004;0.0481000000000003;0.0498999999999983;0.0484999999999998;0.0480999999999983;0.0468999999999995;0.0480999999999999;0.0491000000000005;0.0470999999999991;0.0493000000000006;0.0475999999999989]
 500[5.71360000000000;5.67990000000000;5.71470000000000;5.69170000000000;5.67630000000000;5.67710000000000;5.70790000000000;5.67190000000000;5.68050000000000;5.55540000000000;5.65590000000000;5.60680000000000;5.60390000000000;5.70530000000000;5.63670000000000;5.69100000000000;5.70250000000000;5.68810000000000;5.77050000000000;5.65230000000000;5.70110000000000]
 100[0.0916000000000008;0.0961000000000009;0.0924999999999990;0.0906999999999996;0.108600000000000;0.0961999999999996;0.0912999999999997;0.0901999999999987;0.0809000000000012;0.0905999999999985;0.110900000000000;0.109000000000000;0.0993999999999996;0.108400000000000;0.124499999999999;0.0980999999999996;0.106400000000000;0.105500000000001;0.125200000000001;0.131000000000001;0.116900000000000]
 predict
 1000:[0.0549000000000004;0.0554999999999996;0.0537000000000005;0.0550000000000003;0.0553000000000002;0.0535999999999994;0.0545999999999997;0.0551999999999993;0.0541999999999994;0.0542000000000004;0.0547999999999999;0.0540999999999991;0.0553999999999990;0.0537000000000002;0.0548000000000012;0.0564999999999995;0.0546000000000001;0.0561999999999995;0.0536999999999998;0.0573000000000006;0.0526000000000009]
 500:[6.43740000000000;6.21690000000000;6.30430000000000;6.29950000000000;6.15490000000000;6.16170000000000;6.06240000000000;6.12240000000000;6.07150000000000;6.20780000000000;6.09360000000000;6.10650000000000;6.03810000000000;5.97490000000000;6.03400000000000;6.02800000000000;5.96670000000000;6.02900000000000;6.00610000000000;5.94440000000000;5.86210000000000;0.0549000000000004]
 100:[0.0964;0.0937000000000000;0.0906999999999996;0.0930000000000007;0.111300000000000;0.0953000000000003;0.0941999999999997;0.0891000000000002;0.0819000000000003;0.0914999999999992;0.111500000000000;0.115300000000000;0.0981999999999996;0.112700000000002;0.121699999999998;0.103300000000000;0.106400000000000;0.111700000000001;0.125600000000001;0.129500000000000;0.115000000000000;6.43740000000000]