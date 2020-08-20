addpath('../FM');
addpath('../util');
[A,rows,cols,entries,rep,field,symm]=mmread('../util/well1033.mtx');
b=ones(rows,1);
x0=zeros(cols,1);
maxIter=300;
      [xkR,xkR2,countFR,countNWR,bNWR,tfR,vkR]=residualR(x0,A,b,maxIter);
         dR=norm(xkR-xs);
         rkR=b-A*xkR;
         rkR(rkR<0)=0;
         gR=norm(A'*rkR);
         dR2=norm(xkR2-xs);
         rkR2=b-A*xkR2;
         rkR2(rkR2<0)=0;
         gR2=norm(A'*rkR2);
         fprintf('resdual$ %d \\times %d $ & %g & %g & %4.5f & %d & %d &\n',m,n,gR2,gR,tfR,countFR,countNWR);
         
         [xkA,rkA,countFA,countNA,bNWA,tfA,vkA,Arr]=als(x0,A,b,maxIter);
         dA=norm(xkA-xs);
         rkA=b-A*xkA;
         rkA(rkA<0)=0;
         dA=norm(rkA);
         gA=norm(A'*rkA);
         fprintf('ALS$ %d \\times %d $ & %g & %g & %4.5f & %d & %d & %d\n',m,n,dA,gA,tfA,countFA,countNA,Arr(1,end));
%         Arecord=[Arecord;m n Arr(1,end) gA];
        
%          [xkpa,rkpa,countFMpa,countNWpa,beginNWpa,tfpa,vkpa]=pina(x0,A,b,maxIter);
%          dpa=norm(xkpa-xs);
%          rkpa=b-A*xkpa;
%          rkpa(rkpa<0)=0;
%          gpa=norm(A'*rkpa);
%          fprintf('pina$ %d \\times %d $ & %g & %g & %4.2f & %d & %d &\n',m,n,dpa,gpa,tfpa,countFMpa,countNWpa);
         
         [xkhan,rkhan,countFMhan,countNWhan,beginNWhan,tfhan,vkhan]=han(x0,A,b,maxIter);
         dhan=norm(xkhan-xs);
         rkhan=b-A*xkhan;
         rkhan(rkhan<0)=0;
         dhan=norm(rkhan);
         ghan=norm(A'*rkhan);
         fprintf('han$ %d \\times %d $ & %g & %g & %4.5f & %d & %d &\n',m,n,dhan,ghan,tfhan,countFMhan,countNWhan);
                  
         [xkLei,rkLei,countFMLei,countNWLei,beginNWLei,tfLei,vkLei]=Lei(x0,A,b,maxIter);
         dLei=norm(xkLei-xs);
         rkLei=b-A*xkLei;
         rkLei(rkLei<0)=0;
          dLei=norm(rkLei);
         gLei=norm(A'*rkLei);
         fprintf('Lei$ %d \\times %d $ & %g & %g & %4.5f & %d & %d &\n',m,n,dLei,gLei,tfLei,countFMLei,countNWLei);
         
        [xkD,rkD,countFD,countND,bNWD,tfD,vkD]=Dax(x0,A,b,maxIter);
          xs=-1;
         dD=norm(xkD-xs);
         rkD=b-A*xkD;
         rkD(rkD<0)=0;
         dD=norm(rkD);
         gD=norm(A'*rkD);
         fprintf('Dax$ %d \\times %d $ & %g & %g & %4.5f & %d & %d & %d\n',m,n,dD,gD,tfD,countFD,countND,bNWD);
        
        [xkG,rkG,countFG,countNG,bNWG,tfG,vkG]=gradientFM_i(x0,A,b,nf,1e-8,maxIter,xs,2);
          dG=norm(xkG-xs);
          rkG=b-A*xkG;
          rkG(rkG<0)=0;
          dG=norm(rkG);
          gG=norm(A'*rkG);
          fprintf('grad$ %d \\times %d $ & %g & %g & %4.2f & %g & %g & %g &\n',m,n,dG,gG,tfG,countFG,countNG,bNWG);
 
     [xkC,rkC,countFMC,countNWC,beginNWC,tfC,vkC]=contraction_i(x0,A,b,nf,0.8,maxIter,xs,2);
         dC=norm(xkC-xs);
        rkC=b-A*xkC;
          rkC(rkC<0)=0;
          dC=norm(rkC);
         gC=norm(A'*rkC);
          fprintf('con$ %d \\times %d $ & %g & %g & %4.5f & %g & %g & %g &\n',m,n,dC,gC,tfC,countFMC,countNWC,beginNWC);

       [xkP,rkP,countFP,countNP,bNWP,tfP,vkP]=predictFM_i(x0,A,b,nf,10,maxIter,xs,2);
        dP=norm(xkP-xs);
        rkP=b-A*xkP;
          rkP(rkP<0)=0;
          dP=norm(rkP);
        gP=norm(A'*rkP);
         fprintf('pred$ %d \\times %d $ & %g & %g & %4.5f & %g & %g & %g &\n',m,n,dP,gP,tfP,countFP,countNP,bNWP); 
         record=[m,n,gR,tfR,countFR,countNWR,nf;...
                 m,n,gA,tfA,countFA,countNA,nf;...
                 m,n,gD,tfD,countFD,countND,nf;...
                 m,n,gG,tfG,countFG,countNG,nf;...
                 m,n,gC,tfC,countFMC,countNWC,nf;...
                 m,n,gP,tfP,countFP,countNP,nf...
                ];

        