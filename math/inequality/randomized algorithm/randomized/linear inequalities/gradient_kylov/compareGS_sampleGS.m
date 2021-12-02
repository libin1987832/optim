clear
clc
debug = 1;
%% 产生问题矩阵
% 随机矩阵
m = 100;
n = 20;
r = 50;
m = 700;
n = 2*r+1;

%% generate nodes
%tj are drawing randomly from a uniform distribution in [0,1]
t = rand(m,1);

%then ordering by magnitude (since they are all positive, just sort)
t = sort(t);
% just to assign the size
w = zeros(m,1); 
x = zeros(n,1);
%% generate x 
realx = randn(n,1);
imgx = randn(n,1);
for l = 1:n
   x(l) = realx(l)+1i*imgx(l);
end  
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

% A = 2 * rand(m , n)-1;
b = 2 * rand(m , 1)-1;
% b=A*ones(n,1);
x0 = zeros(n , 1);
% save('test2.mat','A','b','x0')
%load('test.mat');
% 二维矩阵
% A = -[1,-1;-1,-1;0,1];b=-[0;-1;0];x0=[-1;0];
% 不一致情况下的正解
% x_exact=[1/2;1/3];


%% 基于IFM的算法找到一个解
maxit_LSQR = 3;
r = b - A * x0;
r(r<0) = 0;
norm_r0 = norm(r);
norm_g0 = norm(A'*r);
fprintf('%s & %g & %g \n','最开始的目标函数和梯度', norm_r0, norm_g0);
[x_exact, ~, ~, ~, ~] = IFM(A, b, x0,10, maxit_LSQR , 1e-10,[],debug);
r = b - A * x_exact;
r(r<0) = 0;
norm_rexact = norm(r);
norm_gexact = norm(A'*r);
fprintf('%s & %g & %g \n','IFM解的目标函数值和梯度  ', norm_rexact, norm_gexact);
x_exact=[];
%% 参数的设定
 maxit_IFM = 200;
% 
tol=1e-10;
tol=[];


%% GuassSeidel
maxit_Rand =10000;
t=clock;
 [x_GS,iter_GS,error_GS,xA_GS,index_GS] = GuassSeidelNE(A, b, x0,2.0,maxit_Rand,tol,x_exact,debug);
tf_GS=etime(clock,t);
r = b - A * x_GS;
r(r<0) = 0;
r_GS = norm(r);
g_GS = norm(A'*r);
fprintf('& %s & %g & %g & %d & %g \\\\\n', 'GuassSeidel', r_GS, g_GS, iter_GS, tf_GS);

%maxit_Rand =350000;
t=clock;
[x_WGS,iter_WGS,error_WGS,xA_WGS,index_WGS] = randomizedGaussSeidelsampleNE(A, b, x0,2.0,maxit_Rand,tol,x_exact,debug);

tf_WGS=etime(clock,t);
r = b - A * x_WGS;
r(r<0) = 0;
r_WGS = norm(r);
g_WGS = norm(A'*r);
fprintf('& %s & %g & %g & %d & %g \\\\\n', 'weight Gauss', r_WGS, g_WGS, iter_WGS, tf_WGS);
%% 画图
if debug
figure
% h=semilogy(xA_IFM, error_IFM, 'k.');
% h.LineStyle = '--';
display=1:1:maxit_Rand;
h=semilogy(xA_GS(display), error_GS(display), 'r+');
h.LineStyle = '--';
hold on
h=semilogy(xA_WGS(display), error_WGS(display), 'b*');
h.LineStyle = '--';
legend('Gauss Seidel','Rand Guass');
xlabel('the iterative numbers');
ylabel('the norm of the gradient');
end
Acol=sum(A.*A,1);
weight = Acol/sum(Acol);
tabulate(index_GS(2:end))
tt = tabulate(index_WGS(2:end));
[tt(:,3) [weight]'*100]


