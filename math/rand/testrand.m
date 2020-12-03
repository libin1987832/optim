clc
clear
m = 3000;
n = 20;
A = rand(m, n);
x = rand(n, 1);
b = A*x;
x0 = zeros(n, 1);
tol = 1e-13;

maxit = 1000;
AA = A'*A;
tf = issymmetric(AA);
d = eig(AA);
tol1 = length(d)*eps(max(d));
isposdef = all(d > tol1)
%issemidef = all(d > -tol)
[x,fl0,rr0,it0,rv0] = pcg(A'*A,A'*b,tol,maxit);
[x0, rv1, xv1,randp] = randkaczmarz(A, b, x0, tol, maxit);
Axv1 = A * xv1 - repmat(b, 1, size(xv1,2));
xl = xv1 - repmat(x,1,size(xv1,2));
xxs = sum(xl .* xl, 1);

% [d,out]=lineData(A,b,[-1,1],[-1,1]);
% hold on 
% line(d(out,[1,2])',d(out,[3,4])')
% c = num2str([1:10]);
% for i = 1:10
% line([xv1(1,i) xv1(1,i+1)],[xv1(2,i),xv1(2,i+1)])
% text(xv1(1,i), xv1(2,i),num2str(i));
% end

% arrayfun(@(a) text(a.X,a.Y,a.C),s);



semilogy(0:length(rv0)-1,rv0/norm(b),'-o')
hold on
semilogy(0:length(rv1)-1,rv1/norm(b),'-o')
% yline(tol,'r--'); 2018b
legend('pcg','kaczmarz','Tolerance','Location','East')
xlabel('Iteration number')
ylabel('Relative residual')
