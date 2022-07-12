clc
clear

filename='bcsstk19.mtx';
[A,rows,cols,entries,rep,field,symm]=mmread(filename);
[~, n] = size(A)
min(eig(A))
I=eye(n);
B=I;  %初始矩阵B

A = 0.5 * ( A + A' );
B = 0.5 * ( B + B' );

%% 计算初始向量
r=zeros(n,1);
for i=1:n
 r(i)=min(A(:,i)*B(i,i)-A(i,i)*B(:,i));
 if r(i)>=0
  disp('该矩阵不需要')
 end
end 
warning off 
s=find(r==max(r));
m=zeros(n,1);
m(s)=1;
x0=m; %初始向量x1
% x0=unidrnd(2,n,1);
%  x0=rand(n,1);
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
% tic;
% [x, iter] = spBas(A, B, x0, 1e-5, unifrnd (0,1), eps, maxIt);
% toc
% lambda = (x' * A * x) / (x' * B * x);
% disp(['BAS:lambda=' num2str(lambda)  ',Dualfeas='  num2str(min((A - lambda * B) * x ),'%.4e') ',iter=' num2str(iter)])

%% SPL
tic;
%  A1=A+R1*B;
[x, iter]=SPL(A, B, x0, maxIt,  eps, epsbbp);
toc
lambda = (x' * A * x) / (x' * B * x);
disp(['spl:lambda=' num2str(lambda,'%.4e')  ',Dualfeas='  num2str(min((A - lambda * B) * x ),'%.4e') ',iter=' num2str(iter)])
%% FQP
% tic;
% %  A1=A+R1*B;
% [x,iter] = FqpEicp(B, A, x0, maxIt, maxItflqp, maxItsub, 1, eps, epsflqp, epsbbp);
% toc
% lambda = (x' * A * x) / (x' * B * x);
% disp(['FQP_BBP:lambda=' num2str(lambda)  ',Dualfeas='  num2str(min((A - lambda * B) * x ),'%.4e') ',iter=' num2str(iter)])



3

 
%% SSQP(D)
 M=diag(diag(A));
%   M=diag(diag(A+R1*I));%A对称非正定
 
tic;
[x,  iter, error] = sqpEicp(A, B, M, x0, sigma0, 0.1, eps, epsbisect, maxIt, 0);
toc
lambda = (x' * A * x) / (x' * B * x);
disp(['sqp2:lambda=' num2str(lambda,'%.4e')  ',Dualfeas='  num2str(min((A - lambda * B) * x ),'%.4e') ',iter=' num2str(iter)])
%% SQP(G)
M=A; %A对称正定
%    M=A+R1*I;%A对称非正定

tic;
[x,  iter, error] = sqpEicp(A, B, M, x0, sigma0, 0.1, eps, epsbisect, maxIt, 0);
toc
lambda = (x' * A * x) / (x' * B * x);
disp(['sqp1:lambda=' num2str(lambda,'%.4e')  ',Dualfeas='  num2str(min((A - lambda * B) * x ),'%.4e') ',iter=' num2str(iter)])


%% SSQP(I)
M=I;

tic;
[x,  iter, error] = sqpEicp(A, B, M, x0, sigma0, 0.1, eps, epsbisect, maxIt, 0);
toc
lambda = (x' * A * x) / (x' * B * x);
disp(['sqp3:lambda=' num2str(lambda,'%.4e')  ',Dualfeas='  num2str(min((A - lambda * B) * x ),'%.4e') ',iter=' num2str(iter)])

%%
% tic;
% [x,iter] = FqpEicp(B, A, x0, 1000, 10000, 200, 0, 1e-6, 1e-10, 1e-8);
% toc
% lambda = (x' * A * x) / (x' * B * x);
% disp(['FQP_QUAD:lambda=' num2str(lambda) ', ninfx=' num2str(sum(x<epsx)) ',ninfy='  num2str(sum((A - lambda * B) * x < epsxlambda)) ',iter=' num2str(iter)])
% 6

