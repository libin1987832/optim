clear
clc
debug = 1;
%% 产生问题矩阵
m = 8000;
n = 800;
A = 2 * rand(m , n)-1;
b = 2 * rand(m , 1)-1;

x0 = zeros(n , 1);
tol=1e-10-1;
maxit_Rand =500000;
maxit_PC = 50;
maxit_Han = 50;
 maxit_IFM =1000;
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

z0=zeros(m,1);
t=clock;
[x_PC,z_PC,iter_PC,error_PC,index_PC]=PC(x0,z0,A,b,maxit_PC,tol,x_exact,debug);
tf_PC=etime(clock,t);
r = b - A * x_PC;
r(r<0) = 0;
r_PC = norm(r);
g_PC = norm(A'*r);
fprintf('& %s & %g & %g & %d & %g \\\\\n','PC', r_PC, g_PC,iter_PC,tf_PC);

%% 参数的设定

z0=zeros(m,1);
t=clock;
[x_Han,r_Han,countFM_Han,error_Han,beginNW_Han,tf_Han,vk_Han,rkArr_Han]=han(x0,A,b,maxit_Han,1e-11,x_exact,debug);
tf_PC=etime(clock,t);
r = b - A * x_Han;
r(r<0) = 0;
r_PC = norm(r);
g_PC = norm(A'*r);
fprintf('& %s & %g & %g & %d & %g \\\\\n','Han', r_PC, g_PC,countFM_Han,tf_PC);

%% 参数的设定

t=clock;
[x_IFM,iter_IFM,error_IFM,xA_IFM,index_IFM] = IFM(A, b, x0, maxit_IFM, maxit_LSQR ,tol, x_exact,debug);
tf_IFM=etime(clock,t);
r = b - A * x_IFM;
r(r<0) = 0;
r_IFM = norm(r);
g_IFM = norm(A'*r);
fprintf('& %s & %g & %g & %d & %g \\\\\n','IFM', r_IFM, g_IFM,iter_IFM,tf_IFM);

%% GuassSeidel

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
% t=clock;
% [x_GS,iter_RGS,error_RGS,xA_GS,index_GS] = randGuassSeidelNE(A, b, x0,2.0,maxit_Rand*5,tol,x_exact,debug);
% tf_GS=etime(clock,t);
% r = b - A * x_GS;
% r(r<0) = 0;
% r_GS = norm(r);
% g_GS = norm(A'*r);
% fprintf('& %s & %g & %g & %d & %g \\\\\n', 'RCD', r_GS, g_GS, iter_RGS, tf_GS);
% 最开始的目标函数和梯度 & 4.75884 & 9.87503 
% & PC & 4.23416 & 4.41305e-15 & 200000 & 3.844 \\
% & Han & 4.23416 & 2.72853e-15 & 3 & 0.002 \\
% & IFM & 4.23416 & 1.48122e-15 & 10000 & 0.329 \\
% & CCD & 4.23416 & 4.20128e-12 & 100000 & 0.251 \\
% & UCD & 4.23416 & 2.90361e-13 & 100000 & 0.349 \\
% & RCD & 4.23416 & 3.17314e-14 & 500000 & 1.6 \\
% 矩阵维数 & 500 & 50 tol=-1
% 最开始的目标函数和梯度 & 8.36887 & 39.1119 
% & PC & 7.09282 & 7.54872e-13 & 300000 & 13.98 \\
% & Han & 7.09282 & 1.40796e-14 & 4 & 0.009 \\
% & IFM & 7.09282 & 8.44486e-15 & 10000 & 0.813 \\
% & CCD & 7.09282 & 4.59006e-13 & 100000 & 0.327 \\
% & UCD & 7.09282 & 5.51679e-13 & 100000 & 0.444 \\
% & RCD & 7.09282 & 1.18217e-12 & 500000 & 2.164 \\
% 
% 矩阵维数 & 500 & 50 tol=-1
% 最开始的目标函数和梯度 & 8.88904 & 30.9545 
% & PC & 8.22951 & 6.09693e-13 & 500000 & 22.872 \\
% & Han & 8.22951 & 1.93735e-14 & 3 & 0.006 \\
% & IFM & 8.22951 & 1.18153e-14 & 10000 & 0.794 \\
% & CCD & 8.22951 & 2.17866e-12 & 100000 & 0.326 \\
% & UCD & 8.22951 & 1.53127e-12 & 100000 & 0.416 \\
% & RCD & 8.22951 & 5.50179e-12 & 500000 & 2.142 \\
% 
% 矩阵维数 & 500 & 50 tol=-1
% 最开始的目标函数和梯度 & 9.25503 & 34.3172 
% & PC & 8.47067 & 5.51149e-13 & 700000 & 31.928 \\
% & Han & 8.47067 & 1.56159e-14 & 4 & 0.007 \\
% & IFM & 8.47067 & 1.11499e-14 & 10000 & 0.778 \\
% & CCD & 8.47067 & 1.16292e-12 & 100000 & 0.333 \\
% & UCD & 8.47067 & 1.05741e-12 & 100000 & 0.417 \\
% & RCD & 8.47067 & 3.22922e-12 & 500000 & 2.171 \\

% & PC & 8.08725 & 1.94019 & 1000 & 0.054 \\
% & Han & 8.07965 & 1.9053e-14 & 4 & 0.006 \\
% & IFM & 8.07965 & 1.08837e-14 & 30000 & 2.373 \\
% & CCD & 8.07965 & 4.80282e-12 & 200000 & 0.702 \\
% & UCD & 8.07965 & 2.51635e-12 & 200000 & 0.844 \\
% & RCD & 8.07965 & 3.09829e-12 & 1000000 & 4.461 \\

% 矩阵维数 & 1000 & 100 tol=-1
% 最开始的目标函数和梯度 & 12.3856 & 62.4764 
% & PC & 11.2963 & 2.5957e-12 & 200000 & 17.652 \\
% & Han & 11.2963 & 3.83384e-14 & 4 & 0.03 \\
% & IFM & 11.2963 & 2.603e-14 & 30000 & 4.934 \\
% & CCD & 11.2963 & 4.51018e-12 & 200000 & 0.855 \\
% & UCD & 11.2963 & 3.4833e-12 & 200000 & 1.063 \\
% & RCD & 11.2963 & 8.35663e-12 & 1000000 & 5.568 \\

% 矩阵维数 & 1000 & 100 tol=-1
% 最开始的目标函数和梯度 & 13.1045 & 83.4502 
% & PC & 11.4232 & 3.15711e-12 & 500000 & 45.833 \\
% & Han & 11.4232 & 3.78117e-14 & 4 & 0.033 \\
% & IFM & 11.4232 & 2.35302e-14 & 3000 & 0.507 \\
% & CCD & 11.4232 & 3.96318e-12 & 200000 & 0.912 \\
% & UCD & 11.4232 & 2.94397e-12 & 200000 & 1.087 \\
% & RCD & 11.4232 & 6.26025e-12 & 1000000 & 5.846 \\

% & PC & 16.0685 & 1.35482e-11 & 500000 & 161.179 \\
% & Han & 16.0685 & 1.26994e-13 & 4 & 0.915 \\
% & IFM & 16.0685 & 7.62385e-14 & 3000 & 1.782 \\
% & CCD & 16.0685 & 1.31796e-11 & 200000 & 1.48 \\
% & UCD & 16.0685 & 1.03375e-11 & 200000 & 1.766 \\
% & RCD & 16.0685 & 2.79865e-11 & 1000000 & 8.886 \\
% 矩阵维数 & 2000 & 200 tol=-1
% 最开始的目标函数和梯度 & 17.8044 & 140.196 
% & PC & 16.0817 & 1.40915e-11 & 300000 & 92.057 \\
% & Han & 16.0817 & 1.26351e-13 & 4 & 0.801 \\
% & IFM & 16.0817 & 7.55674e-14 & 3000 & 1.798 \\
% & CCD & 16.0817 & 1.6652e-11 & 400000 & 2.966 \\
% & UCD & 16.0817 & 1.54257e-11 & 400000 & 3.396 \\
% 
% 矩阵维数 & 4000 & 400 tol=-1
% 最开始的目标函数和梯度 & 25.5754 & 283.622 
% & PC & 23.6182 & 134.437 & 1000 & 3.591 \\
% & Han & 22.9093 & 3.24553e-13 & 4 & 7.171 \\
% & IFM & 22.9093 & 2.1565e-13 & 3000 & 23.07 \\
% & CCD & 22.9093 & 5.06549e-11 & 400000 & 9.945 \\
% & UCD & 22.9093 & 4.34705e-11 & 400000 & 10.223 \\
% 矩阵维数 & 6000 & 600 tol=-1
% 最开始的目标函数和梯度 & 31.555 & 463.656 
% & PC & 29.6933 & 295.178 & 50 & 0.474 \\
% & Han & 28.0009 & 5.77107e-13 & 4 & 24.633 \\
% & IFM & 28.0009 & 3.95972e-13 & 1000 & 17.885 \\
% & CCD & 28.0009 & 7.71999e-11 & 400000 & 13.603 \\