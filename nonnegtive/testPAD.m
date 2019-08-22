addpath(genpath(pwd));

A1=[1,-1,-1;-2,1,-1;1,1,-1;1,1,-1];
b=[2;-1;-2;3];
n=3; 
m1=4;

iter=3;

fb=b;
fb(fb<0)=0;
%初值
x30=zeros(n,1)+1;
%初值的误差 用来推测优化的结果是不是比这个值还要小
init=0.5*fb'*fb


y=b-A1*x30;
fb(y<0)=0;
f3=0.5*fb'*fb
nec3=[];
All=1:m1;
NE=All(y>0.0001);
All(NE)=0;
nec3=[nec3;NE];
[x31,f3]=PAD(x30,A1,b,NE,iter)
yy=b-A1*x31


