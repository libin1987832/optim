% %°¸Àý1
%  e=0.0001;
%  m1=500;
% m2=500;
% n=700;
% A1=sprand(m1,n,0.1,1/100);
% A2=sprand(m2,n,0.1,1/100);
% b1=rand(m1,1);
% b2=rand(m2,1);
% A=[A1;-A2];
% b=[b1;-b2];
% x0=ones(n,1);
% det=ones(n,1);
% A2(A2<0)=0.5;
% A1(A1<0)=1;
% %°¸Àý2
% A=[1 2 3;4 5 6;-7 -8 -9];
% b=[10 11 -12]';

% x0=[1 1 1]';
% e=0.00001;
% det=[1 1 1]';
% 
% %
addpath('./IPG')
addpath('./exact')
addpath('./hybrid')
clear 
clc
m1 = 200; 
m2 = 200; n = 10;
A1=sprand(m1,n,0.1,1/100);
A2=sprand(m2,n,0.1,1/100);
b1=rand(m1,1);
b2=rand(m2,1);
A=[A1;-A2];
b=[b1;-b2];
x0=ones(n,1);
det=ones(n,1);
maxIterA = 50;
beginp = 1;
options = optimoptions('lsqlin','Algorithm','interior-point','Display','iter');
options.Display = 'off';
options.StepTolerance = 1e-4;
[xk1, resvec, arvec,face1v,face2v, tf]=fixedMatrix(A,b,x0,maxIterA,1e-15,options);
% h=semilogy(beginp:maxIterA,arvec(beginp:maxIterA),'b+');
% h.LineStyle = '--';
% hold on 
alpha = min(eig(A'*A));
[xk2,resvec,arvec,face1h,face2h,tf] = hybridnnls(A,b,x0,alpha,5,maxIterA);
[rpk1, normr1, xmin1, Ar, normKKT1 , face11, face21] = kktResidual(A, b, xk1 , [], 1);
[rpk2, normr2, xmin2, Ar, normKKT2 , face12, face22] = kktResidual(A, b, xk2 , [], 1);
[norm(xk1-xk2)]
[normKKT1 normKKT2;normr1 normr2]

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
