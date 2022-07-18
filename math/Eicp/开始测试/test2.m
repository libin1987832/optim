%%求出当A为对称正定矩阵，B为对称矩阵的特征值互补问题.
%% 随机生成n阶初始矩阵A和B
clc
clear
%% 给定精度
eps=1e-5;
epsbisect=1e-12;
maxIt=100000;

%%
n=4;
k=100;
%  A=[0 1 0 0 0 0;1 0 0 1 1 1;0 0 0 1 0 0;0 1 1 0 1 1; 0 1 0 1 0 0; 0 1 0 1 0 0];
 A=[4 -7 0 0 ;-7 -2 6 0; 0 6 2 -1; 0 0 -1 0];
I=eye(n);
B=I;
%A=A+100*B
L=zeros(k,1);
i=0;
O=inv(B)*A;
R=max(abs(eig(O)));
sigma0 = R + 0.01;
R1=max(abs(eig(A)))+0.01;%A的谱半径  
epsbbp=1e-8;
%%
% x0=[1/6;1/6;1/6;1/6;1/6;1/6];
% [x, iter]=SPL(A, B, x0, maxIt, eps, epsbbp);
% lambda = (x' * A * x) / (x' * B * x)
M=diag(diag(A+R1*I));
%  M=A+R1*I;
% I=eye(n);
% M=I;
 while 1
  i=i+1;
  x0=rand(n,1);
  [x,  iter, error] = sqpEicp(A, B, M, x0, sigma0, 0.1, eps, epsbisect, maxIt, 0);
  lambda = (x' * A * x) / (x' * B * x);
  L(i)= lambda;
  if i==k
     break
  end
 end

L
