clc
clear
addpath('../FM')
addpath('../IFM')
addpath('../strategies')
addpath('../hybrid')
addpath('../linearsolve')

m = 1000;
ratio=0.5;
n = floor( ratio * m);
A = 2 * rand(m , n)-1;
b = 2 * rand(m , 1)-1;
x0 = zeros(n , 1);
maxIter = 200;
nf = 10;
[xkD,rkD,countFD,countND,bNWD,tfD,vkD,rkArrD]=sdax(x0,A,b,maxIter);

rkD=b-A*xkD;
rkD(rkD<0)=0;
dD=norm(rkD);
gD=norm(A'*rkD);
fprintf('Dax$ %d \\times %d $ & %g & %g & %4.2f & %d & %d & %d\n',m,n,dD,gD,tfD,countFD,countND,bNWD);

[xkG,rkG,countFG,countNG,bNWG,tfG,vkG,rkArrG]=sgradientFM(x0,A,b,nf,1e-8,maxIter);
rkG=b-A*xkG;
rkG(rkG<0)=0;
dG=norm(rkG);
gG=norm(A'*rkG);
fprintf('grad$ %d \\times %d $ & %g & %g & %4.2f & %g & %g & %g &\n',m,n,dG,gG,tfG,countFG,countNG,bNWG);

% [xkC,rkC,countFMC,countNWC,beginNWC,tfC,vkC,rkArrC]=scontraction(x0,A,b,nf,0.8,maxIter);
  [xkC,rkC,countFMC,countNWC,beginNWC,tfC,vkC,rkArrC]=contraction_i(x0,A,b,nf,0.8,maxIter,-1,1);
rkC=b-A*xkC;
rkC(rkC<0)=0;
dC=norm(rkC);
gC=norm(A'*rkC);
fprintf('con$ %d \\times %d $ & %g & %g & %4.2f & %g & %g & %g &\n',m,n,dC,gC,tfC,countFMC,countNWC,beginNWC);

[xkP,rkP,countFP,countNP,bNWP,tfP,vkP,rkArrP]=spredictFM(x0,A,b,nf,3,maxIter);

   
rkP=b-A*xkP;
rkP(rkP<0)=0;
dP=norm(rkP);
gP=norm(A'*rkP);
fprintf('pred$ %d \\times %d $ & %g & %g & %4.2f & %g & %g & %g &\n',m,n,dP,gP,tfP,countFP,countNP,bNWP);
