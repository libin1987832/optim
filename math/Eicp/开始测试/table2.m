%%求出当A为对称正定矩阵，B为对称矩阵的特征值互补问题.
%% 随机生成n阶初始矩阵A和B
clc
clear
n=4;
A=[4 -7 0 0 ;-7 -2 6 0; 0 6 2 -1; 0 0 -1 0]
I=eye(n);
B=I;  %初始矩阵B

A = 0.5 * ( A + A' );
B = 0.5 * ( B + B' );

% x1=unidrnd(2,n,1)
x0=[0,0,1,0]';
%% 给定精度
eps=1e-5;
epsbisect=1e-12;
epsbbp=1e-8;
epsflqp=1e-6;
maxIt=100000;
maxItflqp=3000;
maxItsub=3000;

epsx = 0;
epsxlambda = -1e-3;


O=inv(B)*A;
R=max(abs(eig(O)));
sigma0 = R + 0.01;
R1=max(abs(eig(A)))+0.01;%A的谱半径  

%% BAS
tic;
[x, iter] = spBas(A, B, x0, 1e-5, unifrnd (0,1), eps, maxIt);
toc
lambda = (x' * A * x) / (x' * B * x);
disp(['BAS:lambda=' num2str(lambda)  ',Dualfeas='  num2str(min((A - lambda * B) * x ),'%.4e') ',iter=' num2str(iter)])

%% SPL
tic;
A1=A+R1*B;
[x, iter]=SPL(A1, B, x0, maxIt,  eps, epsbbp);
toc
lambda = (x' * A * x) / (x' * B * x);
disp(['spl:lambda=' num2str(lambda)  ',Dualfeas='  num2str(min((A - lambda * B) * x ),'%.4e') ',iter=' num2str(iter)])
%% FQP
tic;
A1=A+R1*B;
[x,iter] = FqpEicp(B, A1, x0, maxIt, maxItflqp, maxItsub, 1, eps, epsflqp, epsbbp);
toc
lambda = (x' * A * x) / (x' * B * x);
disp(['FQP_BBP:lambda=' num2str(lambda)  ',Dualfeas='  num2str(min((A - lambda * B) * x ),'%.4e') ',iter=' num2str(iter)])



3
%% SQP(G)
% M=A; %A对称正定
  M=A+R1*I;%A对称非正定

tic;
[x,  iter, error] = sqpEicp(A, B, M, x0, sigma0, 0.1, eps, epsbisect, maxIt, 0);
toc
lambda = (x' * A * x) / (x' * B * x);
disp(['sqp1:lambda=' num2str(lambda,'%.4e')  ',Dualfeas='  num2str(min((A - lambda * B) * x ),'%.4e') ',iter=' num2str(iter)])


 
%% SSQP(D)
% M=diag(diag(A));
M=diag(diag(A+R1*I));%A对称非正定

tic;
[x,  iter, error] = sqpEicp(A, B, M, x0, sigma0, 0.1, eps, epsbisect, maxIt, 0);
toc
lambda = (x' * A * x) / (x' * B * x);
disp(['sqp2:lambda=' num2str(lambda,'%.4e')  ',Dualfeas='  num2str(min((A - lambda * B) * x ),'%.4e') ',iter=' num2str(iter)])

%% SSQP(I)
M=I;

tic;
[x,  iter, error] = sqpEicp(A, B, M, x0, sigma0, 0.1, eps, epsbisect, maxIt, 0);
toc
lambda = (x' * A * x) / (x' * B * x);
disp(['sqp3:lambda=' num2str(lambda,'%.4e')  ',Dualfeas='  num2str(min((A - lambda * B) * x ),'%.4e') ',iter=' num2str(iter)])

%%
% tic;
% [x,iter] = FqpEicp(B, A, x1, 1000, 10000, 200, 0, 1e-6, 1e-10, 1e-8);
% toc
% lambda = (x' * A * x) / (x' * B * x);
% disp(['FQP_QUAD:lambda=' num2str(lambda) ', ninfx=' num2str(sum(x<epsx)) ',ninfy='  num2str(sum((A - lambda * B) * x < epsxlambda)) ',iter=' num2str(iter)])
% 6

