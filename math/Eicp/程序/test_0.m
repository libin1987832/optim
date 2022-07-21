%%  Test Problem0
clc
clear
n=4;
A=[4 -7 0 0 ;-7 -2 6 0; 0 6 2 -1; 0 0 -1 0]; %The initial matrix A
I=eye(n);
B=I;  %The initial matrix B
A = 0.5 * ( A + A' );
B = 0.5 * ( B + B' );

x0=[0,0,1,0]'; % the initial  point
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
R1=max(abs(eig(A)))+0.01;
%% SQP(G)
M=A+R1*I;

tic;
[x,  iter, error] = sqpEicp(A, B, M, x0, sigma0, 0.1, eps, epsbisect, maxIt, 0);
toc
lambda = (x' * A * x) / (x' * B * x);
disp(['SQP(G):lambda=' num2str(lambda,'%.4e')  ',Dualfeas='  num2str(min((A - lambda * B) * x ),'%.4e') ',iter=' num2str(iter)])
 
%% SSQP(D)
M=diag(diag(A+R1*I));

tic;
[x,  iter, error] = sqpEicp(A, B, M, x0, sigma0, 0.1, eps, epsbisect, maxIt, 0);
toc
lambda = (x' * A * x) / (x' * B * x);
disp(['SSQP(D):lambda=' num2str(lambda,'%.4e')  ',Dualfeas='  num2str(min((A - lambda * B) * x ),'%.4e') ',iter=' num2str(iter)])

%% SSQP(I)
M=I;

tic;
[x,  iter, error] = sqpEicp(A, B, M, x0, sigma0, 0.1, eps, epsbisect, maxIt, 0);
toc
lambda = (x' * A * x) / (x' * B * x);
disp(['SSQP(I):lambda=' num2str(lambda,'%.4e')  ',Dualfeas='  num2str(min((A - lambda * B) * x ),'%.4e') ',iter=' num2str(iter)])

