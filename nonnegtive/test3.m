addpath(genpath(pwd));

A1=[1,-1;-2,1;1,1;1,1];
b=[2;-1;-2;3];
n=2; 
m1=4;
e=0.01;

% % % number test
% m1=10;m2=600;n=3;density=1;cond=10;delt=1e-5;e=1e-4;
% A1=sprand(m1,n,density,1/cond); 
% b=sprand(m1,1,density,1/cond)-0.2;

iter=3;

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
for k=1:41
%记录每次迭代后积极面情况
y=b-A1*x10;
All=1:m1;
NE=All(y>0.0001);
All(NE)=0;
nec1=[nec1;All];
%前面这些步骤主要是用来分析积极面情况
[x11,f1]=FM(x10,Q,R,A1,b);
%记录本次目标函数值
fc1=[fc1,f1];

%记录每次迭代后积极面情况
y=b-A1*x20;
All=1:m1;
NE=All(y>0.0001);
All(NE)=0;
nec2=[nec2;All];
%前面这些步骤主要是用来分析积极面情况
[x21,f2]=IFM(x20,A1,b,iter);
%记录本次目标函数值
fc2=[fc2,f2];

y=b-A1*x30;
All=1:m1;
NE=All(y>0.0001);
All(NE)=0;
nec3=[nec3;All];
[x31,f3]=PAD(x30,A1,b,NE,iter);

fc3=[fc3,f3];
% 跟新迭代值
x10=x11;
x20=x21;
x30=x31;
end
% 局部解
disp "迭代N次后，固定，不精确，积极三种方法最终值"
xn=[x10,x20,x30]
% 局部解
disp "迭代N次后，固定，不精确，积极三种方法最终函数值"
% [norm(b-A1*x10),norm(b-A1*x20),norm(b-A1*x30)]
fcc=[fc1;fc2;fc3]
disp "积极面(0)"
[nec1(end,:);nec2(end,:);nec3(end,:)]
nec3;
fb=b-A1*x30;
fb(fb<0)=0;
A1'*fb
