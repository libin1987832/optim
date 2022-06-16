                                                                                                                                                                                     %%求出当A为对称（正定）矩阵，B为对称正定矩阵的特征值互补问题.
%% 随机生成n阶初始矩阵A和B
clc
clear
n=200;
C=unidrnd(10,n,n);
[C,Y]=qr(C,0);
a=1;
b=1;
z1=a+b.*unidrnd(100,n,1);
Q=diag(z1);
A=C'*Q*C; %初始矩阵A

C=-8+22*rand(n,n);
I=eye(n);
B=C'*C+I;%初始矩阵BC'*C+

%% 计算初始向量
r=zeros(n,1);
for i=1:n
 r(i)=min(A(:,i)*B(i,i)-A(i,i)*B(:,i));
 if r(i)>=0
  disp('该矩阵不需要')
 end
end 
s=find(r==max(r));
m=zeros(n,1);
m(s)=1;
x1=m; %初始向量x1
%% 调用所有函数

% [i]=bas1B(A,x1,B,n); %调用BAS函数，输出时间和迭代次数
% 
 [iG]=SQPGB(A,x1,B,n,I); %调用SQP(G)函数，输出时间和迭代次数
  [x, error] = sqpEicp(A, B, M, x0, sigma0, eps, bisectEps, maxIT, debug)
% 
% [iD]=SSQPDB(A,x1,B,n,I); %调用SSQP(D)函数，输出时间和迭代次数
% 
% [iI]=SSQPIB(A,x1,B,n); %调用SSQP(I)函数，输出时间和迭代次数
% % % % 
% [xk,i,lamdab]=SPL(A,B,x1,n);

