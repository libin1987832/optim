clear
clc
debug = 1;
%% 产生问题矩阵
m = 1000;
n = 200;
% A = 2 * rand(m , n)-1;

r = 100;
m = 6000;
n = 2*r+1;
% generate nodes
%tj are drawing randomly from a uniform distribution in [0,1]
t = rand(m,1);

%then ordering by magnitude (since they are all positive, just sort)
t = sort(t);
%just to assign the size
w = zeros(m,1); 
x = zeros(n,1);
% generate x 
realx = randn(n,1);
imgx = randn(n,1);
for l = 1:n
   x(l) = realx(l)+1i*imgx(l);
end  
% generate A and b
A = zeros(m,n);
for j = 1:m
%    dealing with special cases when reach the endpoints(nodes)
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

A=real(A);
b = 2 * rand(m , 1)-1;
x0 = zeros(n , 1);
tol=1e-10;
fprintf('%s & %d & %d tol=%g\n','矩阵维数', m, n, tol);
%% 基于IFM的算法找到一个解
maxit_LSQR = 3;
r = b - A * x0;
r(r<0) = 0;
norm_r0 = norm(r);
norm_g0 = norm(A'*r);
fprintf('%s & %g & %g \n','最开始的目标函数和梯度', norm_r0, norm_g0);
% [x_exact, ~, ~, ~, ~] = IFM(A, b, x0,10, maxit_LSQR , 1e-10,[],debug);
% r = b - A * x_exact;
% r(r<0) = 0;
% norm_rexact = norm(r);
% norm_gexact = norm(A'*r);
% fprintf('%s & %g & %g \n','IFM解的目标函数值和梯度  ', norm_rexact, norm_gexact);
 x_exact=[];
%% 参数的设定
maxit_PC = 1000;
z0=zeros(m,1);
t=clock;
[x_PC,z_PC,iter_PC,error_PC,index_PC]=PC(x0,z0,A,b,maxit_PC,tol,x_exact,debug);
tf_PC=etime(clock,t);
r = b - A * x_PC;
r(r<0) = 0;
r_PC = norm(r);
g_PC = norm(A'*r);
fprintf('& %s & %g & %g & %d & %g \\\\\n','PC', r_PC, g_PC,iter_PC,tf_PC);

%% GuassSeidel
maxit_Rand =10000;
t=clock;
 [x_GS,iter_GS,error_GS,xA_GS,index_GS] = GuassSeidelNE(A, b, x0,2.0,maxit_Rand,tol,x_exact,debug);
tf_GS=etime(clock,t);
r = b - A * x_GS;
r(r<0) = 0;
r_GS = norm(r);
g_GS = norm(A'*r);
fprintf('& %s & %g & %g & %d & %g \\\\\n', 'CCD', r_GS, g_GS, iter_GS, tf_GS);


%% simpleGuassSeidel
% maxit_Rand =1000000;
t=clock;
[x_GS,iter_SGS,error_SGS,xA_GS,index_GS] = simpleGuassSeidelNE(A, b, x0,2.0,maxit_Rand,tol,x_exact,debug);
tf_GS=etime(clock,t);
r = b - A * x_GS;
r(r<0) = 0;
r_GS = norm(r);
g_GS = norm(A'*r);
fprintf('& %s & %g & %g & %d & %g \\\\\n', 'UCD', r_GS, g_GS, iter_GS, tf_GS);

%% randGuassSeidel
% maxit_Rand =630000;
t=clock;
[x_GS,iter_RGS,error_RGS,xA_GS,index_GS] = randGuassSeidelNE(A, b, x0,2.0,maxit_Rand,tol,x_exact,debug);
tf_GS=etime(clock,t);
r = b - A * x_GS;
r(r<0) = 0;
r_GS = norm(r);
g_GS = norm(A'*r);
fprintf('& %s & %g & %g & %d & %g \\\\\n', 'RCD', r_GS, g_GS, iter_GS, tf_GS);

% 
figure (1)
iters=0;
semilogy(1:iters,error_PC(1:iters),'k--');
title('比较坐标下降算法和投影收缩算法的性能差异')
ylabel('梯度的范数') 
xlabel('迭代次数（对于坐标下降法，每n次迭代算一次迭代）') 
hold on

semilogy(1:iters,error_GS(1:iters),'b--');
semilogy(1:iters,error_SGS(1:iters),'r--');
semilogy(1:iters,error_RGS(1:iters),'g--');

% semilogy(1:size(error_GS,2),error_GS,'b--');
% semilogy(1:size(error_SGS,2),error_SGS,'r--');
% semilogy(1:size(error_RGS,2),error_RGS,'g--');
hold off
legend('PC','CCD','UCD','RCD')


% 矩阵维数 & 1000 & 200 tol=1e-10
% 最开始的目标函数和梯度 & 12.7924 & 108.799 
% & PC & 11.2041 & 69.4224 & 10 & 0.005 \\
% & CCD & 11.4511 & 73.3405 & 100 & 0.004 \\
% & UCD & 11.5413 & 72.8916 & 100 & 0.004 \\
% & RCD & 11.5747 & 73.6909 & 100 & 0.006 \\

% 矩阵维数 & 1000 & 200 tol=1e-10
% 最开始的目标函数和梯度 & 12.7797 & 109.965 
% & PC & 9.76072 & 10.4459 & 2000 & 0.31 \\
% & CCD & 9.68613 & 2.74816 & 1000 & 0.008 \\
% & UCD & 9.75077 & 10.3931 & 1000 & 0.01 \\
% & RCD & 9.80608 & 15.2416 & 1000 & 0.012 \\

% 矩阵维数 & 1000 & 200 tol=1e-10
% 最开始的目标函数和梯度 & 13.1598 & 101.044 
% & PC & 10.5524 & 5.68376 & 3000 & 0.412 \\
% & CCD & 10.5248 & 0.267343 & 2000 & 0.012 \\
% & UCD & 10.5355 & 3.77691 & 2000 & 0.016 \\
% & RCD & 10.5339 & 3.5738 & 2000 & 0.02 \\


% 矩阵维数 & 1000 & 200 tol=1e-10
% 最开始的目标函数和梯度 & 13.3502 & 123.286 
% & PC & 9.77148 & 3.68011 & 4000 & 0.505 \\
% & CCD & 9.75938 & 0.0271003 & 3000 & 0.013 \\
% & UCD & 9.76036 & 1.12242 & 3000 & 0.017 \\
% & RCD & 9.75982 & 0.776219 & 3000 & 0.02 \\

% 
% 矩阵维数 & 1000 & 200 tol=1e-10
% 最开始的目标函数和梯度 & 12.9318 & 97.0274 
% & PC & 10.2512 & 0.161184 & 10000 & 1.152 \\
% & CCD & 10.2512 & 0.00022609 & 5200 & 0.021 \\
% & UCD & 10.2512 & 0.0382373 & 5200 & 0.031 \\
% & RCD & 10.2512 & 0.0520124 & 5200 & 0.034 \\
% y1=[0.005 0.31 0.505 1.152];
% y2=[0.004 0.008 0.013 0.021];
% y3=[0.004 0.01 0.017 0.031];
% y4=[0.006 0.012 0.02 0.034];
% x1=[69.4224,10.4459,3.68011,0.161184];
% x2=[73.3405,2.74816,0.0271003,0.00022609];
% x3=[72.8916,10.3931,1.12242,0.0382373 ];
% x4=[73.6909,15.2416,0.776219,0.0520124];
% x1=[0.005 0.31 0.505 1.152];
% x2=[0.004 0.008 0.013 0.021];
% x3=[0.004 0.01 0.017 0.031];
% x4=[0.006 0.012 0.02 0.034];
% y1=[69.4224,10.4459,3.68011,0.161184];
% y2=[73.3405,2.74816,0.0271003,0.00022609];
% y3=[72.8916,10.3931,1.12242,0.0382373 ];
% y4=[73.6909,15.2416,0.776219,0.0520124];
% figure (1)
% semilogy(x1,y1,'k--');
% title('比较坐标下降算法和投影收缩算法的精度和时间关系')
% ylabel('梯度的范数') 
% xlabel('运行的时间') 
% hold on
% iters=10;
% semilogy(x2,y2,'b--');
% semilogy(x3,y3,'r--');
% semilogy(x4,y4,'g--');
% 
% % semilogy(1:size(error_GS,2),error_GS,'b--');
% % semilogy(1:size(error_SGS,2),error_SGS,'r--');
% % semilogy(1:size(error_RGS,2),error_RGS,'g--');
% hold off
% legend('PC','CCD','UCD','RCD')
% set(gca,'yscale','log')
% set(gca,'xscale','log')