clc
clear
filename='nos4.mtx';
[A,rows,cols,entries,rep,field,symm]=mmread(filename);
[~, n] = size(A);
min(eig(A))
I=eye(n);
B=I;  

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
R1=max(abs(eig(A)))+0.01; 

%% SPL
t0=clock;
[x, iter]=SPL(A, B, x0, maxIt,  eps, epsbbp); %When A is a symmetric positive definite matrix
%A1=A+R1*B; % When A is a symmetric indefinite matrix
%[x, iter]=SPL(A1, B, x0, maxIt,  eps, epsbbp);
tf=etime(clock,t0);
lambda = (x' * A * x) / (x' * B * x);
disp(['spl:lambda=' num2str(lambda,'%.4e')  ',Dualfeas='  num2str(min((A - lambda * B) * x ),'%.4e') ',iter=' num2str(iter) ',CPU=' num2str(tf,'%.4e')])
%% SSQP(D)
M=diag(diag(A)); %When A is a symmetric positive definite matrix
% M=diag(diag(A+R1*I));%When A is a symmetric indefinite matrix
 
t0=clock;
[x,  iter, error] = sqpEicp(A, B, M, x0, sigma0, 0.1, eps, epsbisect, maxIt, 0);
tf=etime(clock,t0);
lambda = (x' * A * x) / (x' * B * x);
disp(['SSQP(D):lambda=' num2str(lambda,'%.4e')  ',Dualfeas='  num2str(min((A - lambda * B) * x ),'%.4e') ',iter=' num2str(iter) ',CPU=' num2str(tf,'%.4e')])
 
%% SQP(G)
M=A;   %When A is a symmetric positive definite matrix
% M=A+R1*I;%When A is a symmetric indefinite matrix

t0=clock;
[x,  iter, error] = sqpEicp(A, B, M, x0, sigma0, 0.1, eps, epsbisect, maxIt, 0);
tf=etime(clock,t0);
lambda = (x' * A * x) / (x' * B * x);
disp(['SQP(G):lambda=' num2str(lambda,'%.4e')  ',Dualfeas='  num2str(min((A - lambda * B) * x ),'%.4e') ',iter=' num2str(iter) ',CPU=' num2str(tf,'%.4e')])
 
%% SSQP(I)
M=I;

t0=clock;
[x,  iter, error] = sqpEicp(A, B, M, x0, sigma0, 0.1, eps, epsbisect, maxIt, 0);
tf=etime(clock,t0);
lambda = (x' * A * x) / (x' * B * x);
disp(['SSQP(I):lambda=' num2str(lambda,'%.4e')  ',Dualfeas='  num2str(min((A - lambda * B) * x ),'%.4e') ',iter=' num2str(iter) ',CPU=' num2str(tf,'%.4e')])
 