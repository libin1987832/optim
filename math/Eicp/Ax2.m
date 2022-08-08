%%求出当A为对称（正定）矩阵，B为对称正定矩阵的特征值互补问题.
%% 随机生成n阶初始矩阵A和B
clc
clear

%% 给定初始矩阵A和B
n=500;
C=unidrnd(10,n,n);
[C,Y]=qr(C,0);
a=100;
b=1;
z1=0:(a-b)/(n-1):1;
z1=z1+b;
for i=1:n
    z(i)=b+((i-1)*(a-b))/(n-1); %在区间（a，b）内产生均匀分布的n维向量
end
Q=diag(z);
A=C'*Q*C; %初始矩阵A
A=(A+A')/2;
C=-8+22*rand(n,n);
I=eye(n);
B=I;%初始矩阵BC'*C+
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
% x1=[1; 1; 1 ; 1;0];
%% 给出罚参数
O=inv(B)*A;
R=max(abs(eig(O))); 
sigma0 = R + 0.01;

%% 调用
R1=max(abs(eig(A)))+0.01;

 tic 
 [iD]=SSQPDB(A,x1,B,n,I); %调用SSQP(D)函数，输出时间和迭代次数
toc

% M=A; %A对称正定 
% % M=A+R1*I; %A对称非正定 
% typeM = 1;
% method = 1 ;
% tic;
% [x,  iter] = sqpEicp(A, B, M, x1, sigma0, 0.1, 1e-5, 1e-12, 10000, typeM, method );
% toc
% lambda = (x' * A * x) / (x' * B * x); 
% disp(['lambda=' num2str(lambda) ', ninfx=' num2str(sum(x<0)) ',ninfy='  num2str(sum((A - lambda * B) * x < -1e-3)) ',iter=' num2str(iter)])

M=diag(diag(A));%A对称正定
% M=diag(diag(A+R1*I)); %A对称非正定 
typeM =3; 
method =2 ;
tic;
[x,  iter ] = sqpEicp(A, B, M, x1, sigma0, 0.1, 1e-5, 1e-12, 10000, typeM, method);
toc
lambda = (x' * A * x) / (x' * B * x);
disp(['lambda=' num2str(lambda) ', ninfx=' num2str(sum(x<0)) ',ninfy='  num2str(sum((A - lambda * B) * x < -1e-3)) ',iter=' num2str(iter)])


% M=I;
% typeM = 3;
% method = 2;
% tic;
% [x,  iter] = sqpEicp(A, B, M, x1, sigma0, 0.1, 1e-5, 1e-12, 10000, typeM, method);
% toc
% lambda = (x' * A * x) / (x' * B * x);
% disp(['lambda=' num2str(lambda) ', ninfx=' num2str(sum(x<0)) ',ninfy='  num2str(sum((A - lambda * B) * x < -1e-3)) ',iter=' num2str(iter)])


tic;
[x, iter, error]=SPL(A, B, x1, 10000,  1e-5, 1e-6, 0);
toc
lambda = (x' * A * x) / (x' * B * x);
disp(['lambda=' num2str(lambda) ', ninfx=' num2str(sum(x<0)) ',ninfy='  num2str(sum((A - lambda * B) * x < -1e-3)) ',iter=' num2str(iter)])


% tic;
% [x, iter, fun] = spBas(A, B, x1, 1e-5, unifrnd (0,1), 1e-5, 10000, 0);
% toc
% lambda = (x' * A * x) / (x' * B * x);
% disp(['lambda=' num2str(lambda) ', ninfx=' num2str(sum(x<0)) ',ninfy='  num2str(sum((A - lambda * B) * x < -1e-3)) ',iter=' num2str(iter)])
% 
% [i]=bas1B(A,x1,B,n); %调用BAS函数，输出时间和迭代次数
