addpath('../FM');
addpath('../IFM');
addpath('../util');
addpath('../strategies');
addpath('../linearsolve');
% [A,rows,cols,entries,rep,field,symm]=mmread('../util/well1033.mtx');
[A,rows,cols,entries,rep,field,symm]=mmread('../util/well1850.mtx');
m=rows;
n=cols;
% A(20:20:end,:)=0;
b=rand(rows,1);
% b(1:2:end)=-1;
% m=1000;
% n=100;
%  A=2*rand(m,n)-1;
%  b=2*rand(m,1)-1;
 x0=ones(n,1)*100000;
maxIter=1000;
nf=2;
type=1;
r0=b-A*x0;
r0(r0<0)=0;
dr0=norm(A'*r0);
fprintf('resdual$ dim       $ & ||f(x)|| & ||df(x)|| & time & iteration & begin subspace %g&\n',dr0);
      [xkR,xkR2,countFR,countNWR,bNWR,tfR,vkR]=residualR(x0,A,b,maxIter);
      xs=xkR;
         dR=norm(xkR-xs);
         rkR=b-A*xkR;
         rkR(rkR<0)=0;
         gR=norm(A'*rkR);
          dR2=norm(xkR2-xs);
         rkR2=b-A*xkR2;
         rkR2(rkR2<0)=0;
         gR2=norm(A'*rkR2);
         fprintf('resdual$ %d \\times %d $ & %g & %g & %4.5f & %d & %d & %g,%g&\n',m,n,gR2,gR,tfR,countFR,countNWR,sum(rkR>1e-13),sum(rkR>=0));
         
         [xkA,rkA,countFA,countNA,bNWA,tfA,vkA,Arr]=als(x0,A,b,maxIter);
         dA=norm(xkA-xs);
         rkA=b-A*xkA;
         rkA(rkA<0)=0;
          dA=norm(rkA);
         gA=norm(A'*rkA);
         fprintf('ALS$ %d \\times %d $ & %g & %g & %4.5f & %d & %d & %d& %g,%g\n',m,n,dA,gA,tfA,countFA,countNA,Arr(1,end),sum(rkA>1e-13),sum(rkA>=0));
%         Arr(1,end)
%         Arecord=[Arecord;m n Arr(1,end) gA];
        
%          [xkpa,rkpa,countFMpa,countNWpa,beginNWpa,tfpa,vkpa]=pina(x0,A,b,maxIter);
%          dpa=norm(xkpa-xs);
%          rkpa=b-A*xkpa;
%          rkpa(rkpa<0)=0;
%          gpa=norm(A'*rkpa);
%          fprintf('pina$ %d \\times %d $ & %g & %g & %4.2f & %d & %d &\n',m,n,dpa,gpa,tfpa,countFMpa,countNWpa);
         xs1=-1;
         [xkhan,rkhan,countFMhan,countNWhan,beginNWhan,tfhan,vkhan]=han(x0,A,b,maxIter);
         rkhan=b-A*xkhan;
         rkhan(rkhan<0)=0;
         dhan=norm(rkhan);
         ghan=norm(A'*rkhan);
         fprintf('han$ %d \\times %d $ & %g & %g & %4.5f & %d & %d & %g,%g&\n',m,n,dhan,ghan,tfhan,countFMhan,countNWhan,sum(rkhan>1e-13),sum(rkhan>=0));
                  
         [xkLei,rkLei,countFMLei,countNWLei,beginNWLei,tfLei,vkLei]=Lei(x0,A,b,maxIter);
         rkLei=b-A*xkLei;
         rkLei(rkLei<0)=0;
          dLei=norm(rkLei);
         gLei=norm(A'*rkLei);
         fprintf('Lei$ %d \\times %d $ & %g & %g & %4.5f & %d & %d & %g,%g&\n',m,n,dLei,gLei,tfLei,countFMLei,countNWLei,sum(rkLei>1e-13),sum(rkLei>=0));
         
        [xkD,rkD,countFD,countND,bNWD,tfD,vkD]=Dax(x0,A,b,maxIter);
         rkD=b-A*xkD;
         rkD(rkD<0)=0;
          dD=norm(rkD);
         gD=norm(A'*rkD);
         fprintf('Dax$ %d \\times %d $ & %g & %g & %4.5f & %d & %d & %d& %g,%g\n',m,n,dD,gD,tfD,countFD,countND,bNWD,sum(rkD>1e-13),sum(rkD>=0));
        
        [xkG,rkG,countFG,countNG,bNWG,tfG,vkG]=gradientFM_i(x0,A,b,nf,1e-10,maxIter,xs1,type);
          rkG=b-A*xkG;
          rkG(rkG<0)=0;
          dG=norm(rkG);
          gG=norm(A'*rkG);
          fprintf('grad$ %d \\times %d $ & %g & %g & %4.2f & %g & %g & %g & %g,%g&\n',m,n,dG,gG,tfG,countFG,countNG,bNWG,sum(rkG>1e-13),sum(rkG>=0));
 
     [xkC,rkC,countFMC,countNWC,beginNWC,tfC,vkC]=contraction_i(x0,A,b,nf,0.999,maxIter,xs1,type);
        rkC=b-A*xkC;
          rkC(rkC<0)=0;
          dC=norm(rkC);
         gC=norm(A'*rkC);
          fprintf('con$ %d \\times %d $ & %g & %g & %4.5f & %g & %g & %g & %g,%g&\n',m,n,dC,gC,tfC,countFMC,countNWC,beginNWC,sum(rkC>1e-13),sum(rkC>=0));

       [xkP,rkP,countFP,countNP,bNWP,tfP,vkP]=predictFM_i(x0,A,b,nf,10,maxIter,xs1,type);
        rkP=b-A*xkP;
          rkP(rkP<0)=0;
           dP=norm(rkP);
        gP=norm(A'*rkP);
         fprintf('pred$ %d \\times %d $ & %g & %g & %4.5f & %g & %g & %g & %g,%g&\n',m,n,dP,gP,tfP,countFP,countNP,bNWP,sum(rkP>1e-13),sum(rkP>=0)); 
         record=[m,n,gR,tfR,countFR,countNWR,nf;...
                 m,n,gA,tfA,countFA,countNA,nf;...
                 m,n,gD,tfD,countFD,countND,nf;...
                 m,n,gG,tfG,countFG,countNG,nf;...
                 m,n,gC,tfC,countFMC,countNWC,nf;...
                 m,n,gP,tfP,countFP,countNP,nf...
                ];

     fprintf('distence from han:han(%g),resdual(%g),Lei(%g),Dax(%g),reduction(%g),contract(%g),predict(%g)\n',norm(xkhan),norm(xkhan-xkR),norm(xkhan-xkLei),norm(xkhan-xkD),norm(xkhan-xkG),norm(xkhan-xkC),norm(xkhan-xkP))