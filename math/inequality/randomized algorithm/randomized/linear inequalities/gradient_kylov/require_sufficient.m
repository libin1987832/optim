clear
clc
debug = 1;
%% 产生问题矩阵
% 随机矩阵
repeat = 100;
statis = zeros(repeat,1);
for i = 1:repeat
m = 200;
n = 50;

A = 2 * rand(m , n)-1;
b = 2 * rand(m , 1)-1;
% b=A*ones(n,1);
x0 = zeros(n , 1);
% save('test2.mat','A','b','x0')
% load('test2.mat');
% 二维矩阵
% A = -[1,-1;-1,-1;0,1];b=-[0;-1;0];x0=[-1;0];
% 不一致情况下的正解
% x_exact=[1/2;1/3];
maxit_Rand = 2000;
% 
tol=1e-10;
tol=[];
x_exact=[];
t=clock;
[x_WGS,iter_WGS,error_WGS,xA_WGS,index_WGS] = optimGaussSeidel_sufficient(A, b, x0,2.0,maxit_Rand,tol,x_exact,debug);

tf_WGS=etime(clock,t);
r = b - A * x_WGS;
r(r<0) = 0;
r_WGS = norm(r);
g_WGS = norm(A'*r);
fprintf('& %s & %g & %g & %d & %g \\\\\n', 'weight Gauss', r_WGS, g_WGS, iter_WGS, tf_WGS);

sufficinet = index_WGS(2:end,:);
sums=sum(sufficinet,1);
statis(i)=find(sums,1,'last');
end
ts=tabulate(statis);
% ts(:,1:2)
bar(ts(:,2))
%tabulate(statis)