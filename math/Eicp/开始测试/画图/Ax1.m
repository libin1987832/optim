%%求出当A为对称正定矩阵，B为对称矩阵的特征值互补问题.
%% 随机生成n阶初始矩阵A和B
clc
clear
n=100;
A=importdata('A2.txt');
I=eye(n);
 B=I; %初始矩阵B
%B=importdata('B2.txt');
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
[x, OUT1] = spBas(A, B, x0, 1e-5, unifrnd (0,1), eps, maxIt);
toc
lambda = (x' * A * x) / (x' * B * x);

%% SPL
tic;
 A1=A+R1*B;
[x, OUT2]=SPL(A1, B, x0, maxIt,  eps, epsbbp);
toc
lambda = (x' * A * x) / (x' * B * x);

%% FQP
tic;
 A1=A+R1*B;
[x,OUT3] = FqpEicp(B, A1, x0, maxIt, maxItflqp, maxItsub, 1, eps, epsflqp, epsbbp);
toc
lambda = (x' * A * x) / (x' * B * x);



3
%% SQP(G)
% M=A; %A对称正定
M=A+R1*I;%A对称非正定

tic;
[x,  OUT4] = sqpEicp(A, B, M, x0, sigma0, 0.1, eps, epsbisect, maxIt, 0);
toc
lambda = (x' * A * x) / (x' * B * x);


 
%% SSQP(D)
% M=diag(diag(A));
 M=diag(diag(A+R1*I));%A对称非正定

tic;
[x,  OUT5] = sqpEicp(A, B, M, x0, sigma0, 0.1, eps, epsbisect, maxIt, 0);
toc
lambda = (x' * A * x) / (x' * B * x);

%% SSQP(I)
M=I;

tic;
[x,  OUT6] = sqpEicp(A, B, M, x0, sigma0, 0.1, eps, epsbisect, maxIt, 0);
toc
lambda = (x' * A * x) / (x' * B * x);

%% 画图
figure
  semilogy(OUT1,'LineWidth',1.5);
  hold on
 semilogy(OUT2,'LineWidth',1.5);
  hold on
 semilogy(OUT3,'LineWidth',1.5);
 semilogy(OUT4,'LineWidth',1.5);
 semilogy(OUT5,'LineWidth',1.5);
 semilogy(OUT6,'LineWidth',1.5);
% 
 legend('BAS','SPL','FLQP','SQP(G)','SSQP(D)','SSQP(I)')
title('B=I')
axis([0,600,10E-7,10E+2])
xlabel('iterations')
ylabel('norm(dk)')
