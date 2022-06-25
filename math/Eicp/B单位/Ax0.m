%%求出当A为对称（正定）矩阵，B为单位矩阵的特征值互补问题.
%% 随机生成n阶初始矩阵A和B
clc
clear
addpath('../')
n=50;
C=unidrnd(10,n,n);
[C,Y]=qr(C,0);
a=-21;
b=10;
r=a+b.*unidrnd(100,n,1);% 在区间（a，b）内产生均匀分布的n维向量
Q=diag(r);
A=C'*Q*C; %初始矩阵A
B=eye(n);%初始矩阵B
%% 计算初始向量
r=zeros(n,1);
for i=1:n
 r(i)=min(A(:,i)*B(i,i)-A(i,i)*B(:,i));
 if r(i)>=0
  disp('该矩阵不需要')
 end
end
s=find(r==max(r));
s=s(1);
m=zeros(n,1);
m(s)=1;
x1=m; %初始向量x1
%% 调用所有函数


[xk,i,h,lamdab]=SPL(A,B,x1,n);
[iD]=SSQPD(A,x1,n); %调用SSQP(D)函数，输出时间和迭代次数
R=max(abs(eig(A)));    %A的谱半径                                          %用A的1范数来近似A的谱（代价会不会太大了）
e=zeros(n,1);
e(1)=R+0.01;
I=eye(n);
M=A;%A对称正定
M=A+(R+0.01)*I;
M = diag(diag(M));
 sigma0 = 0;
epsx = 0;
 epsxlambda = -1e-3;
% tic;[x,  iter, error] = sqpEicp(A, B, M, x1, sigma0, 0.1, 1e-5, 1e-12, 10000, 0);toc
% lambda = (x' * A * x) / (x' * B * x);
% disp(['sqp:lambda=' num2str(lambda) ', ninfx=' num2str(sum(x<epsx)) ',ninfy='  num2str(sum((A - lambda * B) * x < epsxlambda)) ',iter=' num2str(iter)])
% [iI]=SSQPI(A,x1,n); %调用SSQP(I)函数，输出时间和迭代次数
% tic;[x,  iter, error] = sqpEicp(A, B, I, x1, sigma0, 0.1, 1e-5, 1e-12, 10000, 0);toc
% lambda = (x' * A * x) / (x' * B * x);
% disp(['sqp:lambda=' num2str(lambda) ', ninfx=' num2str(sum(x<epsx)) ',ninfy='  num2str(sum((A - lambda * B) * x < epsxlambda)) ',iter=' num2str(iter)])
[iG]=SQPG(A,x1,n); %调用SQP(G)函数，输出时间和迭代次数
M=A+(R+0.01)*I ;  %A对称非正定
M=(M+M')/2;
tic;[x,  iter, error] = sqpEicp(A, B, M, x1, sigma0, 0.1, 1e-5, 1e-12, 10000, 0);toc
lambda = (x' * A * x) / (x' * B * x);
disp(['sqp:lambda=' num2str(lambda) ', ninfx=' num2str(sum(x<epsx)) ',ninfy='  num2str(sum((A - lambda * B) * x < epsxlambda)) ',iter=' num2str(iter)])
[i]=bas1(A,x1,n); %调用BAS函数，输出时间和迭代次数
