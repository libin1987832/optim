
%%  Test Problem 2——4
clc
clear
n = input('Performance of algorithms for solving Test Problems, input n='); % Order of the matrix
C=unidrnd(10,n,n);
[C,Y]=qr(C,0);
%% The initial matrix A （the symmetric positive definite matrix）
a=1;
b=1;
z1=a+b.*unidrnd(1000,n,1);
Q=diag(z1);
A=C'*Q*C; 
%% The initial matrix A （the symmetric indefinite definite matrix）
% a=-100;
% b=1;
% z1=a+b.*unidrnd(1000,n,1);
% a1=-unidrnd(20,1,1);
% z1(1)=a1;
% Q=diag(z1);
% A=C'*Q*C; %初始矩阵A
% min(eig(A))
%%  The initial matrix B
C=-8+22*rand(n,n);
I=eye(n);
B=I;         %Identity matrix
%B=I+C'*C;   %The symmetric positive definite matrix

A = 0.5 * ( A + A' );
B = 0.5 * ( B + B' );

%% the initial  point
r=zeros(n,1);
for i=1:n
 r(i)=min(A(:,i)*B(i,i)-A(i,i)*B(:,i));
 if r(i)>=0
  disp('This matrix does not require')  %该矩阵不需要
 end
end 
warning off 
s=find(r==max(r));
m=zeros(n,1);
m(s)=1;
x0=m; 
%% Tolerance
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
R1=max(abs(eig(A)))+0.01; %A的谱半径  

%% BAS
t0=clock;
[x, iter] = spBas(A, B, x0, 1e-5, unifrnd (0,1), eps, maxIt);
tf=etime(clock,t0);
lambda = (x' * A * x) / (x' * B * x);
disp(['BAS:lambda=' num2str(lambda)  ',Dualfeas='  num2str(min((A - lambda * B) * x ),'%.4e') ',iter=' num2str(iter) ',CPU=' num2str(tf,'%.4e')])

%% SPL
t0=clock;
[x, iter]=SPL(A, B, x0, maxIt,  eps, epsbbp); %When A is a symmetric positive definite matrix
%A1=A+R1*B; % When A is a symmetric indefinite matrix
%[x, iter]=SPL(A1, B, x0, maxIt,  eps, epsbbp);
tf=etime(clock,t0);
lambda = (x' * A * x) / (x' * B * x);
disp(['spl:lambda=' num2str(lambda)  ',Dualfeas='  num2str(min((A - lambda * B) * x ),'%.4e') ',iter=' num2str(iter) ',CPU=' num2str(tf,'%.4e')])
%% FLQP
t0=clock;
[x,iter] = FlqpEicp(B, A, x0, maxIt, maxItflqp, maxItsub, 1, eps, epsflqp, epsbbp); %When A is a symmetric positive definite matrix
% A1=A+R1*B; %When A is a symmetric indefinite matrix
% [x,iter] = FqpEicp(B, A1, x0, maxIt, maxItflqp, maxItsub, 1, eps, epsflqp, epsbbp);
tf=etime(clock,t0);
lambda = (x' * A * x) / (x' * B * x);
disp(['FLQP_BBP:lambda=' num2str(lambda)  ',Dualfeas='  num2str(min((A - lambda * B) * x ),'%.4e') ',iter=' num2str(iter) ',CPU=' num2str(tf,'%.4e')])
%% SQP(G)

M=A; %When A is a symmetric positive definite matrix
% M=A+R1*I;%When A is a symmetric indefinite matrix

t0=clock;
[x,  iter, error] = sqpEicp(A, B, M, x0, sigma0, 0.1, eps, epsbisect, maxIt, 0);
tf=etime(clock,t0);
lambda = (x' * A * x) / (x' * B * x);
disp(['SQP(G):lambda=' num2str(lambda)  ',Dualfeas='  num2str(min((A - lambda * B) * x ),'%.4e') ',iter=' num2str(iter) ',CPU=' num2str(tf,'%.4e')])


 
%% SSQP(D)
 M=diag(diag(A));  %When A is a symmetric positive definite matrix
% M=diag(diag(A+R1*I));%When A is a symmetric indefinite matrix
 
t0=clock;
[x,  iter, error] = sqpEicp(A, B, M, x0, sigma0, 0.1, eps, epsbisect, maxIt, 0);
tf=etime(clock,t0);
lambda = (x' * A * x) / (x' * B * x);
disp(['SSQP(D):lambda=' num2str(lambda)  ',Dualfeas='  num2str(min((A - lambda * B) * x ),'%.4e') ',iter=' num2str(iter) ',CPU=' num2str(tf,'%.4e')])


%% SSQP(I)
M=I;

t0=clock;
[x,  iter, error] = sqpEicp(A, B, M, x0, sigma0, 0.1, eps, epsbisect, maxIt, 0);
tf=etime(clock,t0);
lambda = (x' * A * x) / (x' * B * x);
disp(['SSQP(I):lambda=' num2str(lambda)  ',Dualfeas='  num2str(min((A - lambda * B) * x ),'%.4e') ',iter=' num2str(iter) ',CPU=' num2str(tf,'%.4e')])

