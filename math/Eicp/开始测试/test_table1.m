%%求出当A为对称正定矩阵，B为对称矩阵的特征值互补问题.
%% 随机生成n阶初始矩阵A和B
clc
clear
%n=20;
n = input('Performance of algorithms for solving Test Problems, input n=');
typeB = input('Performance of algorithms for solving Test Problems, input Btype(1:I,other:C^TC+I)=');
C=unidrnd(10,n,n);
[C,Y]=qr(C,0);
%% 矩阵A的特征值随机生成
a=1;
b=1;
z1=a+b.*unidrnd(1000,n,1);
Q=diag(z1);
A=C'*Q*C; %初始矩阵A
%% 
C=-8+22*rand(n,n);
I=eye(n);
if typeB == 1
B=I;
else
B=I+C'*C;  %初始矩阵B
end
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

%% Fmincon

fun = @(x)  (x' * A * x)/ 2;
ceq=@(x) 0.5*x' * B * x-1;
nonlcon=@(x)deal([],ceq(x));
disp(['problem:n=' num2str(n)])

t0=clock;
opts = optimoptions('fmincon','Display','none','Algorithm','sqp','MaxIterations',20000);
x = fmincon(fun, x0, [],[], [], [], zeros(n, 1), [], nonlcon, opts);
tf=etime(clock,t0);
lambda = (x' * A * x) / (x' * B * x);
disp(['fmincon:lambda=' num2str(lambda,'%.4e')  ',Dualfeas='  num2str(min((A - lambda * B) * x ),'%.4e') ',CPU=' num2str(tf,'%.4e')])


%% SQP(G)
M=A; %A对称正定
%  M=A+R1*I;%A对称非正定

t0=clock;
[x,  iter, error] = sqpEicp(A, B, M, x0, sigma0, 0.1, eps, epsbisect, maxIt, 0);
tf=etime(clock,t0);
lambda = (x' * A * x) / (x' * B * x);
disp(['SQP(G):lambda=' num2str(lambda,'%.4e')  ',Dualfeas='  num2str(min((A - lambda * B) * x ),'%.4e') ',iter=' num2str(iter) ',CPU=' num2str(tf,'%.4e')])


 
%% SSQP(D)
M=diag(diag(A));
%  M=diag(diag(A+R1*I));%A对称非正定

t0=clock;
[x,  iter, error] = sqpEicp(A, B, M, x0, sigma0, 0.1, eps, epsbisect, maxIt, 0);
tf=etime(clock,t0);
lambda = (x' * A * x) / (x' * B * x);
disp(['SSQP(D):lambda=' num2str(lambda,'%.4e')  ',Dualfeas='  num2str(min((A - lambda * B) * x ),'%.4e') ',iter=' num2str(iter) ',CPU=' num2str(tf,'%.4e')])

%% SSQP(I)
M=I;

t0=clock;
[x,  iter, error] = sqpEicp(A, B, M, x0, sigma0, 0.1, eps, epsbisect, maxIt, 0);
tf=etime(clock,t0);
lambda = (x' * A * x) / (x' * B * x);
disp(['SSQP(I):lambda=' num2str(lambda,'%.4e')  ',Dualfeas='  num2str(min((A - lambda * B) * x ),'%.4e') ',iter=' num2str(iter) ',CPU=' num2str(tf,'%.4e')])


