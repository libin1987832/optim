addpath(genpath(pwd));
% A=[2,1;1,2];
% b=[5;6];

% ezplot('2*x+y-5',[0,10])
% hold on 
% ezplot('x+2*y-6',[0,10])
% [x,fk]=fsearchx(A,b,[10;1],0.000001,0.0001);
% A=-1*A;
% b=-1*b;
% % 精确算法
% [x,f0]=alg1(-1*A,-1*b,[0;10],0.0001);
% % 不精确算法
%  [x,fk]=inexact(-1*b,-1*A,[10;1],0.5,0.6)
% % 罚函数方法
% M=1000000;% 惩罚因子
% delt=0.00001;% delt数值计算的误差
% e=0.0001; %种植条件 梯度来判断
% x=GNP([10;1],M,delt,e,A,b)
% % [mm,vv]=minValue(x,0.1,100,A,b,@fq,@constrain);
% % mm

% % number test
m1=500;m2=500;n=700;density=0.1;cond=100;delt=1e-5;e=1e-4;
A1=sprand(m1,n,density,1/cond);
A2=sprand(m2,n,density,1/cond);
b1=rand(m1,1);b2=rand(m2,1);
A=[A1;-A2];
b=[b1;-b2];
fq(A,b,zeros(n,1))
delt=0.00001;% delt数值计算的误差
e=0.0001; %种植条件 梯度来判断
 [x,fk]=fsearchx(A,b,ones(n,1)*100,e,delt);
% M=1000000;% 惩罚因子

% [x,fk1]=GNP(ones(n,1)*100,M,delt,e,A,b);
% fk
% fk1