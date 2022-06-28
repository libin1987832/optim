%%求出当A为对称（正定）矩阵，B为对称正定矩阵的特征值互补问题.
%% 随机生成n阶初始矩阵A和B
clc
clear
n=500;
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

 n = 4;
% n = 3;
I=eye(n);
 A=[5 7 6 5;7 10 8 7;6 8 10 9;5 7 9 10];
 B = I;
A = 0.5 * ( A + A' );
B = 0.5 * ( B + B' );
% A = A(2:4,2:4);

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
x1=[0.46;0.63;0.61;0]; %初始向量x1
x1=unidrnd(2,n,1);
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
epsxlambda = -1e-3;
[x, iter, error] = allsqp(A, B, M, x1, sigma0, 0.1, 1e-5, 1e-12, 10000, 0);
x