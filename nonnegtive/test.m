addpath(genpath(pwd));

A=[2,1;-3,-1];
b=[5;-3];
n=2;
e=0.01;
 x=pda(A,b,2);
% [x,f]=inexact(-1*b,-1*A,ones(n,1)*100,0.2,0.1,e);
x
% % % number test
% m1=600;m2=600;n=1200;density=0.1;cond=100;delt=1e-5;e=1e-4;
% A1=sprand(m1,n,density,1/cond);
% A2=sprand(m2,n,density,1/cond);
% b1=rand(m1,1);b2=rand(m2,1);
% A=[A1;-A2];
% b=[b1;-b2];
% % 输出零向量对应的函数值
% fqf(b,A,zeros(n,1))
% delt=0.00001;% delt防止病态参数
% e=0.01; %终止条件 梯度的范数小于这个值就终止
% 
% % 精确搜索法
% %   [x0,f0]=alg1(A,b,ones(n,1)*100,e);
% 
% % 非精确搜索法
% % inexact(b,A,ones(n,1)*100,0.2,0.1,e); 
% 
%  M=100;% 惩罚因子
% % 罚函数法
%   [x,fk1]=GNP(ones(n,1)*100,M,delt,e,-1*A,-1*b);
