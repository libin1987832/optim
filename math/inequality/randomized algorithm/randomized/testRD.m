clear
clc
debug = 0;
%% 产生问题矩阵
% 随机矩阵
m = 1000;
n = 200;
maxit=n;
A = 2 * rand(m , n)-1;
b = 2 * rand(m , 1)-1;
B=A'*A;
Acol = diag(B);
p=2;
At_r = A'*b;

[x,rr,index]=Rand_Gauss_Seidel_R(A, b, At_r,maxit,B,Acol,p,1+n/m);
[x,rg]=Gass_seidel_D(A, b, maxit,Acol,1+n/m);

figure
h=semilogy(1:(maxit+1), rg, 'r+');
h.LineStyle = '--';
hold on
h=semilogy(1:(maxit+1), rr, 'k.');
h.LineStyle = '--';


% h=semilogy(xA_In, error_In, 'b*');

legend('Gauss Seidel','Weight rand');
xlabel('the iterative numbers');
ylabel('the norm of the gradient');