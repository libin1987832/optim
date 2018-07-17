addpath(genpath(pwd));

A=[2,1;-3,-1];
b=[5;-3];
n=2;
e=0.01;

% [x,f]=inexact(b,A,ones(n,1)*100,0.2,0.1,e);

% % % number test
m1=6;m2=600;n=3;density=1;cond=10;delt=1e-5;e=1e-4;
A1=sprand(m1,n,density,1/cond); 
b=sprand(m1,1,density,1/cond)-0.1;
fb=b;
fb(fb<0)=0;
x0=zeros(n,1)+0.1;
init=0.5*fb'*fb
[Q,R]=qr(A1);
[x1,f1]=FM(x0,Q,R,A1,b);
f1

[x1,f1]=IFM(x0,A1,b,2);
f1

y=b-A1*x0;
All=1:m1;
NE=All(y>=0)
[x1,f1]=PAD(x0,A1,b,NE,2);
f1

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
