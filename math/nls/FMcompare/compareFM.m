clear
clc
debug = 0;
%% �����������
% �������

r = 500;
m = 5000;
n = 2*r+1;
% t=n;
% n=m;
% m=t;
%% generate nodes
%tj are drawing randomly from a uniform distribution in [0,1]
t = rand(m,1);

%then ordering by magnitude (since they are all positive, just sort)
t = sort(t);
% just to assign the size
w = zeros(m,1); 

%% generate A and b
A = zeros(m,n);
for j = 1:m
    % dealing with special cases when reach the endpoints(nodes)
    if j == 1
        w(j) = (t(2)-t(end)-1)/2;
    elseif j == m
            w(j) = (t(1)+1-t(m-1))/2;
        else
            w(j) = (t(j+1)-t(j-1))/2;
    end              
    for k = -r:r
       A(j,k+r+1) =sqrt(w(j))*exp(2*pi*1i*k*t(j));
    end
end  

x = zeros(n,1);
%% generate x 
realx = randn(n,1);
imgx = randn(n,1);
for l = 1:n
   x(l) = realx(l)+1i*imgx(l);
end  
b = A * x;
x0 = zeros(n , 1);



x_exact=[];
%% �������趨
 maxit_FM = 1000;
maxit_gs = 1;
tol=1e-1;
% tol=[];
alpha = 1;
t=clock;
[x_FM,iter_FM,error_k,iter_k,index_k] = DFM(A, b, x0, maxit_FM,alpha,maxit_gs,tol, x_exact,debug);
tf_FM = etime(clock,t);
r = b - A * x_FM;
r(r<0) = 0;
r_FM = norm(r);
g_FM = norm(A'*r);
fprintf('& %s & %s & %s & %s & %s \\\\\n', 'alg', 'norm(r_+)', 'norm(Ar_+)', 'iteration', 'time');
fprintf('& %s & %g & %g & %d & %g \\\\\n','FM', r_FM, g_FM, iter_FM*(m+n)*maxit_gs, tf_FM);


