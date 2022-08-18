%%求出当A为对称（正定）矩阵，B为对称正定矩阵的特征值互补问题.
%% 随机生成n阶初始矩阵A和B
clc
clear
n=10;
C=unidrnd(10,n,n);
%C = load('C');
%[C,Y]=qr(C.C,0);
[C,Y]=qr(C,0);
a=100;
b=1;
z1=0:(a-b)/(n-1):1;
z1=z1+b;
for i=1:n
    z(i)=b+((i-1)*(a-b))/(n-1);
end
I=eye(n);
Q=diag(z);% 在区间（a，b）内产生均匀分布的n维向量
A=C'*Q*C; %初始矩阵A
B=I+C'*C; %importdata('B3.txt'); %初始矩阵B

%  n = 4;
 n = 6;
I=eye(n);
%  A=[5 7 6 5;7 10 8 7;6 8 10 9;5 7 9 10];
A =[8 -1 0 0; -1 6 -3 0; 0 -3 4 -2; 0 0 -2 5];
%A = [0 1 0 0 0 0;1 0 1 1 0 0;0 1 0 1 1 0;0 1 1 0 1 0;0 0 1 1 0 1;0 0 0 0 1 0] + 5*eye(n) ;
 B = I;
M = I;
%  M=A; %A对称正定
epsx = 0;
epsxlambda = -1e-3;
x1 = ones(n,1)/n;
O=inv(B)*A;
R=max(abs(eig(O))); %A的谱半径  
sigma0 = R + 0.01;
[xa, iter, error] = allsqp(A, B, M, x1, sigma0, 0.1, 1e-10, 1e-12, 1000, 0);
nx = size(xa,2);
xd = []
for i = 1 : nx
    x = xa(:,i) ./ sum(xa(:,i));
lambda = (x' * A * x) / (x' * B * x)-5;
%xd = [xd [lambda; x]];
xd = [xd lambda];
%disp(['FQP_QUAD:lambda=' num2str(lambda) ', ninfx=' num2str(sum(x<epsx)) ',ninfy='  num2str(sum((A - lambda * B) * x < epsxlambda))])
end
 sort(uniquetol(xd,0.01),'descend')'