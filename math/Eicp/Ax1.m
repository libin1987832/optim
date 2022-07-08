%%求出当A为对称（正定）矩阵，B为对称正定矩阵的特征值互补问题.
%% 随机生成n阶初始矩阵A和B
clc
clear

%C = load('C');
%[C,Y]=qr(C.C,0);

%A=C'*Q*C; %初始矩阵A
%B=I+C'*C; %importdata('B3.txt'); %初始矩阵B
[A,rows,cols,entries,rep,field,symm] = mmread('bcsstk02.mtx');
n=cols;
I=speye(n);
B=I;
A = 0.5 * ( A + A' );
B = 0.5 * ( B + B' );
% A = A(2:4,2:4);


warning off 
%x1=[0.46;0.63;0.61;0]; %初始向量x1
x1=sparse(ones(n,1)./(n+1));
%% 调用所有函数
% [iG]=SQPGB(A,x1,B,n,I); %调用SQP(G)函数，输出时间和迭代次数
% M=A;%A对称正定
% O=inv(B)*A;
% R=max(abs(eig(O))); %A的谱半径  
% sigma0 = R + 0.01;
% [x, error] = sqpEicp(A, B, M, x1, sigma0, 1e-5, 1e-5, 10000, 0);
%[iD]=SSQPDB(A,x1,B,n,I); %调用SSQP(D)函数，输出时间和迭代次数
M=diag(A);%A对称正定
O=inv(B)*A;
R=max(abs(eig(O))); %A的谱半径  
sigma0 = 10 + 0.01;
R1=max(abs(eig(A)))+0.01;
% M=A+R1*I;%A对称非正定
M = zeros(n,n);
M = I;
M = 0.5 * (M+M');
%  M=A; %A对称正定
M=diag(diag(M));
epsx = 0;
epsxlambda = -1e-7;
% [x, iter, error] = allsqp(A, B, M, x1, sigma0, 0.1, 1e-5, 1e-12, 10000, 0);

%tic;[x,  iter, error] = sqpEicp(A, B, M, x1, sigma0, 0.1, 1e-5, 1e-12, 10000, 0);toc
%lambda = (x' * A * x) / (x' * B * x);
%disp(['sqp:lambda=' num2str(lambda) ', ninfx=' num2str(sum(x<epsx)) ',ninfy='  num2str(sum((A - lambda * B) * x < epsxlambda)) ',iter=' num2str(iter)])

tic;[x,eta, iter, error]=SPL(A, B, x1, 1,  1e-6, 1e-16, 0);toc
lambda = (x' * A * x) / (x' * x);
w = A * x - lambda * B * x - eta * sparse(ones(n,1)); 
disp(['spl:lambda=' num2str(lambda) ',dualfeasible='  num2str(min(w)) ',iter=' num2str(iter)])

tic;[x, iter, fun] = spBas(A, B, x1, 1e-8, unifrnd (0,1), 1e-8, 1000, 0);toc
lambda = (x' * A * x) / (x' * B * x);
disp(['BAS:lambda=' num2str(lambda) ', ninfx=' num2str(sum(x<epsx)) ',ninfy='  num2str(sum((A - lambda * B) * x < epsxlambda)) ',iter=' num2str(iter)])
x1 = rand(n ,1);
x1 = x1 ./ sum(x1);
tic;[x,iter] = FqpEicp(A, B, x1, 3, 1000, 200, 0, 1e-16, 1e-10, 1e-8);toc
lambda = (x' * A * x) / (x' * B * x);
disp(['FQP_QUAD:lambda=' num2str(lambda) ', ninfx=' num2str(sum(x<epsx)) ',ninfy='  num2str(sum((A - lambda * B) * x < epsxlambda)) ',iter=' num2str(iter)])
tic;[x,iter] = FqpEicp(A, B, x1, 20,1000, 200, 1, 1e-16, 1e-10, 1e-8);toc
lambda = (x' * A * x) / (x' * B * x);
disp(['FQP_BBP:lambda=' num2str(lambda) ', ninfx=' num2str(sum(x<epsx)) ',ninfy='  num2str(sum((A - lambda * B) * x < epsxlambda)) ',iter=' num2str(iter)])
% [i]=bas1B(A,x1,B,n, unifrnd (0,1)); %调用BAS函数，输出时间和迭代次数
% [iI]=SSQPIB(A,x1,B,n); %调用SSQP(I)函数，输出时间和迭代次数
% [xk,i,h,lamdab]=SPL(A,B,x1,n);