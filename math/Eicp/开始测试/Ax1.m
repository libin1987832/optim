%%求出当A为对称正定矩阵，B为对称矩阵的特征值互补问题.
%% 随机生成n阶初始矩阵A和B
clc
clear
n=6;
C=unidrnd(10,n,n);
[C,Y]=qr(C,0);
%% 随机生成对称非正定矩阵A
a=1;
b=1;
z1=a+b.*unidrnd(1000,n,1);
Q=diag(z1);
A=C'*Q*C; %初始矩阵A
%% 随机生成对称非正定矩阵A
% a=-100;
% b=1;
% z1=a+b.*unidrnd(1000,n,1);
% a1=-unidrnd(20,1,1);
% z1(1)=a1;
% Q=diag(z1);
% A=C'*Q*C; %初始矩阵A
% min(eig(A))
%% 生成对称正定矩阵A 
% a=1000;
% b=1;
% z=zeros(n,1);
% for i=1:n
%     z(i)=b+((i-1)*(a-b))/(n-1);
% end
% Q=diag(z);% 在区间（a，b）内产生均匀分布的n维向量
% A=C'*Q*C; %初始矩阵A
%% 
C=-8+22*rand(n,n);
I=eye(n);
B=I+C'*C;  %初始矩阵B

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
 x1=unidrnd(2,n,1)
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
% M=diag(diag(A+R1*I));%A对称非正定
% 
% tic;
% [x,  iter, error] = sqpEicp(A, B, M, x0, sigma0, 0.1, eps, epsbisect, maxIt, 0);
% toc
% lambda = (x' * A * x) / (x' * B * x);
% disp(['sqp2:lambda=' num2str(lambda,'%.4e')  ',Dualfeas='  num2str(min((A - lambda * B) * x ),'%.4e') ',iter=' num2str(iter)])

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

