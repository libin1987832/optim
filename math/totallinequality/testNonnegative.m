% %°¸Àý1
% %
addpath('./IPG')
addpath('./exact')
addpath('./hybrid')
clear 
clc
m1 = 500; 
m2 = 500; n = 70;
 A1=sprand(m1,n,0.1,1/100);
 A2=sprand(m2,n,0.1,1/100);
%A1=rand(m1,n) + 1 ;
%A2=rand(m2,n) + 1;

% b1=rand(m1,1)+1;
% b2=rand(m2,1)+1;
% A=[A1;-A2];
% b=[b1;-b2];
% x0=ones(n,1);

 load('test2.mat')
 n=700;
%det=ones(n,1);
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
% hold on k1
 alpha = 1/max(eig(A'*A));
 maxit = 200;
 
%[x1,f1,residual,exitflag,output,ff]=lsqlin([A,-eye(m1+m2)],b,[],[],[],[],zeros(n+m1+m2,1),Inf*ones(n+m1+m2,1),x0,options);
    
 %xk0 = x1(1:n);
%[rpk0, normr0, xmin0, Ar0, normKKT0 , face120, face220] = kktResidual(A, b, xk0 , [], 1); 
%[xk1, resvec, arvec, face1vec, face2vec, tf1]=fixedMatrix(A,b,x0,maxit,1e-15,options);
%[rpk1, normr1, xmin1, Ar1, normKKT1 , face12, face22] = kktResidual(A, b, xk1 , [], 1); 
% [xk2,resvec,arvec,face1h,face2h,tf2] = hybridnnls(A,b,x0,2,maxIterA,options);
%[rpk2, normr2, xmin2, Ar2, normKKT2 , face12, face22] = kktResidual(A, b, xk2 , [], 1);
 [xk3, resvec3, arvec3, faceXvec3, tf3]  =IPG(A,b,x0,1e-4,1e-5,0.5,5000);
[rpk3, normr3, xmin3, Ar3, normKKT3 , face13, face23] = kktResidual(A, b, xk3 , [], 1);
%[tf2 normKKT2 normr2
 %tf3 normKKT3 normr3]

%
%fprintf('& %s & %g & %g & %g & %g &%g \\\\\n',['FM'],normr1,xmin1,normKKT1,Ar1,tf1);
  fprintf('& %s & %g & %g & %g & %g &%g \\\\\n',['Hybrid'],normr2,xmin2,normKKT2,Ar2,tf2);
  fprintf('& %s & %g & %g & %g & %g &%g \\\\\n',['IPG'],normr3,xmin3,normKKT2,Ar3,tf3);
      
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
