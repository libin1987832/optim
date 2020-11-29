% %°¸Àý1
% %
addpath('./IPG')
addpath('./exact')
addpath('./hybrid')
clear 
clc
m1 = 500; 
m2 = 500; n = 700;
A1=sprand(m1,n,0.1,1/100);
A2=sprand(m2,n,0.1,1/100);
% A1=rand(m1,n);
% A2=rand(m2,n);

b1=rand(m1,1);
b2=rand(m2,1);
A=[A1;-A2];
b=[b1;-b2];
x0=ones(n,1);
% load('test2.mat')
% n=700;
det=ones(n,1);
 maxIterA = 20;
% beginp = 1;
options = optimoptions('lsqlin','Algorithm','interior-point','Display','iter');
options.Display = 'off';
options.StepTolerance = 1e-13;
% [xk1, resvec, arvec,face1v,face2v, tf1]=fixedMatrix(A,b,x0,maxIterA*5,1e-15,options);
% r = b - A * xk1;
% r(r<0) = 0;
% Ar = A' * r;
% KKT = xk1 .* Ar;

% [rpk1, normr1, xmin1, Ar, normKKT1 , face11, face21] = kktResidual(A, b, xk1 , [], 1);
% [normKKT1 max(abs(KKT)) min(xk1) min(Ar)]
% h=semilogy(beginp:maxIterA,arvec(beginp:maxIterA),'b+');
% h.LineStyle = '--';
% hold on 
 alpha = 1/max(eig(A'*A));
[xk2,resvec,arvec,face1h,face2h,tf2] = hybridnnls(A,b,x0,alpha,2,maxIterA,options);
[rpk2, normr2, xmin2, Ar, normKKT2 , face12, face22] = kktResidual(A, b, xk2 , [], 1);
 [xk3,resvec3,arvec3,face1vec3,face2vec3,tf3]=fsearchx(A,b,x0,1e-5,1e-8,40);
[rpk3, normr3, xmin3, Ar3, normKKT3 , face13, face23] = kktResidual(A, b, xk3 , [], 1);

% 
%  fprintf('& %s & %g & %g & %g & %g &%g \\\\\n',['Hybrid'],normr2,xmin2,normKKT2,Ar,tf2);
      
%[norm(xk1-xk2)]
%[normKKT1 normKKT2;normr1 normr2]
%tf1
% face1 = face1(face1>0);
% face2 = face2(face2>0);
% face1f=fliplr(face1);
% face2f=fliplr(face2);
% face1f(1:10)
% face2f(1:10)

% A2(A2<0)=0.5;
% A1(A1<0)=1;
% A=[2,1;2,-1];
% b=[5;-3];
% [x,fk]=fsearchx(A,b,[1;1],0.000001,0.0001)
