clear
clc
debug = 0;
%% 产生问题矩阵
% 随机矩阵

r = 800;
m = 10000;
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


maxit_GS =70000;
x_exact=[];
%% 参数的设定
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

t=clock;
[x_GS,iter_GS,error_GS,xA_GS,index_GS] = GuassSeidelNE(A, b, x0,1.0,maxit_GS,tol,x_exact,debug);
tf_GS=etime(clock,t);
r = b - A * x_GS;
r(r<0) = 0;
r_GS = norm(r);
g_GS = norm(A'*r);
fprintf('& %s & %g & %g & %d & %g \\\\\n', 'GuassSeidel', r_GS, g_GS, iter_GS, tf_GS);

t=clock;
[x_IFMp,iter_IFMp,error_IFMp,xA_IFMp,index_IFMp] = IFMs(A, b, x0, maxit_IFM, maxit_LSQR ,tol, x_exact,debug);
tf_IFMp=etime(clock,t);
r = b - A * x_IFMp;
r(r<0) = 0;
r_IFMp = norm(r);
g_IFMp = norm(A'*r);
fprintf('& %s & %g & %g & %d & %g \\\\\n','IFM', r_IFMp, g_IFMp,iter_IFMp*(m+n)*maxit_LSQR,tf_IFMp);

t=clock;
 [x_GS,iter_GS,error_GS,xA_GS,index_GS] = randGSNE(A, b, x0,1.0,maxit_GS,tol,x_exact,debug);
tf_GS=etime(clock,t);
r = b - A * x_GS;
r(r<0) = 0;
r_GS = norm(r);
g_GS = norm(A'*r);
fprintf('& %s & %g & %g & %d & %g \\\\\n', 'rand', r_GS, g_GS, iter_GS, tf_GS);

maxit_GS =70000;
t=clock;
 [x_GS,iter_GS,error_GS,xA_GS,index_GS] = urandGSNE(A, b, x0,1.0,maxit_GS,tol,x_exact,debug);
tf_GS=etime(clock,t);
r = b - A * x_GS;
r(r<0) = 0;
r_GS = norm(r);
g_GS = norm(A'*r);
fprintf('& %s & %g & %g & %d & %g \\\\\n', 'urand', r_GS, g_GS, iter_GS, tf_GS);


