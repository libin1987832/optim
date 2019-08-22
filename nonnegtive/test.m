addpath(genpath(pwd));

A=[2,1;-3,-1];
b=[5;-3];
n=2;
e=0.01;

% [x,f]=inexact(b,A,ones(n,1)*100,0.2,0.1,e);

% % % number test
m1=10;m2=600;n=3;density=1;cond=10;delt=1e-5;e=1e-4;
A1=sprand(m1,n,density,1/cond); 
b=sprand(m1,1,density,1/cond)-0.2;
%load matlab.mat
iter=3;
% load matlab_test1.mat
% A1=[1 1;-1 -1];
% b=[1;2]; 
% m1=2;
% y=b-A1*[-sqrt(2)/4;-sqrt(2)/4];
% A1'*y;
% y(y<0)=0;
% gg=0.5*y'*y;

fb=b;
fb(fb<0)=0;
%初值
x0=zeros(n,1)+1;
%初值的误差 用来推测优化的结果是不是比这个值还要小
init=0.5*fb'*fb
%Q R 分解 用于固定迭代算法过程中 不用重复分解
[Q,R]=qr(A1'*A1);
%记录 三种方法的的初值
x10=x0;
x20=x0;
x30=x0;
%记录每次迭代后的目标函数值
fc1=[];
fc2=[];
fc3=[];
%记录每次迭代后积极面情况
nec1=[];
nec2=[];
nec3=[];
for k=1:20
%记录每次迭代后积极面情况
y=b-A1*x10;
All=1:m1;
NE=All(y>=0);
All(NE)=0;
nec1=[nec1;All];
%前面这些步骤主要是用来分析积极面情况
[x11,f1]=FM(x10,Q,R,A1,b);
%记录本次目标函数值
fc1=[fc1,f1];

%记录每次迭代后积极面情况
y=b-A1*x20;
All=1:m1;
NE=All(y>=0);
All(NE)=0;
nec2=[nec2;All];
%前面这些步骤主要是用来分析积极面情况
[x21,f2]=IFM(x20,A1,b,iter);
%记录本次目标函数值
fc2=[fc2,f2];



y=b-A1*x30;
All=1:m1;
NE=All(y>=0);
All(NE)=0;
nec3=[nec3;All];
[x31,f3]=PAD(x30,A1,b,NE,iter);
% y=b-A1*x1;
% A1'*y
% y(y<0)=0;
% gg2=0.5*y'*y
% NE=All(y>=0)
fc3=[fc3,f3];
% 跟新迭代值
x10=x11;
x20=x21;
x30=x31;
end
% 局部解
x10
x20
x30
% 局部解
norm(b-A1*x10)
norm(b-A1*x20)
norm(b-A1*x30)
fcc=[fc1;fc2;fc3]
nec1(end,:)
nec2(end,:)
nec3(end,:)
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
