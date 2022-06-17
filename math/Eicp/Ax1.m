%%求出当A为对称（正定）矩阵，B为对称正定矩阵的特征值互补问题.
%% 随机生成n阶初始矩阵A和B
clc
clear
n=100;
C=importdata('C2.txt');
a=100;
b=1;
z1=0:(a-b)/(n-1):1;
z1=z1+b;
for i=1:n
    z(i)=b+((i-1)*(a-b))/(n-1);
end
Q=diag(z);% 在区间（a，b）内产生均匀分布的n维向量
A=C'*Q*C; %初始矩阵A

I=eye(n);
B=importdata('B3.txt'); %初始矩阵B
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
% x1=unidrnd(2,n,1)
%% 调用所有函数
% [iG]=SQPGB(A,x1,B,n,I); %调用SQP(G)函数，输出时间和迭代次数
% M=A;%A对称正定
% O=inv(B)*A;
% R=max(abs(eig(O))); %A的谱半径  
% sigma0 = R + 0.01;
% [x, error] = sqpEicp(A, B, M, x1, sigma0, 1e-5, 1e-5, 10000, 0);
[iD]=SSQPDB(A,x1,B,n,I); %调用SSQP(D)函数，输出时间和迭代次数
M=diag(A);%A对称正定
O=inv(B)*A;
R=max(abs(eig(O))); %A的谱半径  
sigma0 = R + 0.01;
R1=max(abs(eig(A)))+0.01;
 M=A+R1*I;%A对称非正定
%  M=A; %A对称正定
M=diag(diag(M));
[x, error] = sqpEicp(A, B, M, x1, sigma0, 0.1, 1e-5, 1e-12, 10000, 0);
% [i]=bas1B(A,x1,B,n); %调用BAS函数，输出时间和迭代次数
% [iI]=SSQPIB(A,x1,B,n); %调用SSQP(I)函数，输出时间和迭代次数
% [xk,i,h,lamdab]=SPL(A,B,x1,n);