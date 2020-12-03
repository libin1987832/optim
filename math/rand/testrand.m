m = 10;
n = 2;
A = rand(m, n);
b = rand(m, 1);
x0 = zeros(n, 1);
tol = 1e-13;
tol1 = length(d)*eps(max(d));
maxit = 1000;
AA = A'*A;
tf = issymmetric(AA);
d = eig(AA);
isposdef = all(d > tol1)
%issemidef = all(d > -tol)
[x,fl0,rr0,it0,rv0] = pcg(A'*A,A'*b,tol,maxit);
[x0, rv1, xv1,randp] = randkaczmarz(A, b, x0, tol, maxit);
xl = xv1 - repmat(x,1,size(xv1,2));
xxs = sum(xl .* xl, 1);
semilogy(0:length(rv0)-1,rv0/norm(b),'-o')
hold on
semilogy(0:length(rv1)-1,rv1/norm(b),'-o')
% yline(tol,'r--'); 2018b
legend('pcg','kaczmarz','Tolerance','Location','East')
xlabel('Iteration number')
ylabel('Relative residual')
