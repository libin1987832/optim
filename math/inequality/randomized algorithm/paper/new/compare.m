clear
clc
debug = 0;
%% 产生问题矩阵
% 随机矩阵
% m = 500;
% n = 200;


m = 40000;
n = 6000;
  A = 2 * rand(m , n)-1;
   b = 2 * rand(m , 1)-1;

x0 = zeros(n , 1);

tol=1e-10;

fprintf('%s & %d & %d \n','矩阵维数', m, n);
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
 maxit_IFM = 1000;
t=clock;
[x_IFM,iter_IFM,error_IFM,xA_IFM,index_IFM] = IFM(A, b, x0, maxit_IFM, maxit_LSQR ,tol, x_exact,debug);
tf_IFM=etime(clock,t);
r = b - A * x_IFM;
r(r<0) = 0;
r_IFM = norm(r);
g_IFM = norm(A'*r);
fprintf('& %s & %g & %g & %d & %g \\\\\n','IFM', r_IFM, g_IFM,iter_IFM,tf_IFM);

%% GuassSeidel
maxit_Rand =5000000;
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
[x_GS,iter_GS,error_SGS,xA_GS,index_GS] = simpleGuassSeidelNE(A, b, x0,2.0,maxit_Rand,tol,x_exact,debug);
tf_GS=etime(clock,t);
r = b - A * x_GS;
r(r<0) = 0;
r_GS = norm(r);
g_GS = norm(A'*r);
fprintf('& %s & %g & %g & %d & %g \\\\\n', 'UCD', r_GS, g_GS, iter_GS, tf_GS);

%% randGuassSeidel
% maxit_Rand =630000;
t=clock;
[x_GS,iter_GS,error_RGS,xA_GS,index_GS] = randGuassSeidelNE(A, b, x0,2.0,maxit_Rand,tol,x_exact,debug);
tf_GS=etime(clock,t);
r = b - A * x_GS;
r(r<0) = 0;
r_GS = norm(r);
g_GS = norm(A'*r);
fprintf('& %s & %g & %g & %d & %g \\\\\n', 'RCD', r_GS, g_GS, iter_GS, tf_GS);








%% FM
% maxit =100;

% alpha=1;
% maxit_gs=n;
% t=clock;
% [x_FM,iter_FM,error_k,iter_dFM,index_k] = DFM(A, b, x0, maxit_IFM,alpha,maxit_gs,tol/50, x_exact,debug);
% tf_FM=etime(clock,t);
% r = b - A * x_FM;
% r(r<0) = 0;
% r_FM = norm(r);
% g_FM = norm(A'*r);
% fprintf('& %s & %g & %g & %d & %g \\\\\n','FM', r_FM, g_FM,iter_FM,tf_FM);

%% RFM
% maxit =18;
% tol=[];
% alpha=1;
% maxit_R=n;
% t=clock;
% [x_FM,iter_FM,error_k,iter_dFM,index_k] = RFM(A, b, x0, maxit,alpha,maxit_R, tol, x_exact,debug);
% tf_FM=etime(clock,t);
% r = b - A * x_FM;
% r(r<0) = 0;
% r_FM = norm(r);
% g_FM = norm(A'*r);
% fprintf('& %s & %g & %g & %d & %g \\\\\n','RFM', r_FM, g_FM,iter_FM,tf_FM);


% %% han
% maxIter = 5;
% t=clock;
% [x_GS,rkh,countFMh,countNWh,beginNWh,tfh,vkh,rkArrh]=han(x0,A,b,maxIter);
% tf_GS=etime(clock,t);
% r = b - A * x_GS;
% r(r<0) = 0;
% r_GS = norm(r);
% g_GS = norm(A'*r);
% fprintf('& %s & %g & %g & %d & %g \\\\\n', 'han', r_GS, g_GS, countFMh, tf_GS);


% semilogy(error_GS, 'r') % Matrix A
% xlabel('Iterations') 
% ylabel('Error') 
% hold on 
% semilogy(error_SGS, 'k') % Matrix M
% semilogy(error_RGS, 'g') % Matrix M
% semilogy(error_IFM, 'b') % Matrix M
% legend('CCD', 'UCD', 'RCD','IFM' )

%%%%%%%%%%%%
% 矩阵维数 & 8000 & 1401 
% 最开始的目标函数和梯度 & 36.9619 & 28.2101 
% IFM解的目标函数值和梯度   & 7.83054 & 1.8207 
% & GuassSeidel & 6.85029 & 2.84068 & 43431 & 4.646 \\
% & simpleGuassSeidel & 4.65163 & 0.883448 & 29421 & 3.216 \\
% & randGuassSeidel & 4.64761 & 0.883935 & 25218 & 2.868 \\
% & IFM & 4.02986 & 0.623185 & 29 & 5.687 \\
% & FM & 9.11728 & 4.12671 & 33 & 4.011 \\

%%%%%%%%%%%
% 矩阵维数 & 8000 & 5000 
% 最开始的目标函数和梯度 & 36.7672 & 1502.49 
% & IFM & 9.73235e-07 & 9.15091e-06 & 422 & 76.021 \\
% & CCD & 9.85611e-07 & 8.23254e-06 & 1405000 & 39.619 \\
% & UCD & 9.89213e-07 & 8.46126e-06 & 2265000 & 68.175 \\
% & RCD & 9.97063e-07 & 8.54955e-06 & 2235000 & 93.787 \\

%%%%%%%%%
% 矩阵维数 & 7000 & 5000 
% 最开始的目标函数和梯度 & 33.9312 & 1372.63 
% & IFM & 9.89736e-07 & 1.31695e-05 & 183 & 30.229 \\
% & CCD & 9.14101e-07 & 1.01399e-05 & 725000 & 19.681 \\
% & UCD & 9.37602e-07 & 1.07014e-05 & 1075000 & 30.864 \\
% & RCD & 9.87354e-07 & 1.13674e-05 & 1065000 & 42.464 \\

%%%%%%%%%%%%
% 矩阵维数 & 8000 & 6000 
% 最开始的目标函数和梯度 & 35.6655 & 1576.08 
% & IFM & 9.60808e-07 & 1.54778e-05 & 141 & 31.869 \\
% & CCD & 9.47603e-07 & 1.28784e-05 & 672000 & 19.171 \\
% & UCD & 9.43271e-07 & 1.31044e-05 & 1002000 & 29.871 \\
% & RCD & 9.52272e-07 & 1.33195e-05 & 990000 & 42.817 \\

%%%%%%%%%
% 矩阵维数 & 10000 & 8000 
% 最开始的目标函数和梯度 & 41.131 & 2145.7 
% & IFM & 8.9914e-07 & 1.94051e-05 & 102 & 39.228 \\
% & CCD & 8.77339e-07 & 1.53744e-05 & 680000 & 24.703 \\
% & UCD & 9.1509e-07 & 1.65742e-05 & 992000 & 36.899 \\
% & RCD & 9.82543e-07 & 1.76539e-05 & 976000 & 52.616 \\

%%%%%%%%%%%%%%%
% 矩阵维数 & 12000 & 10000 
% 最开始的目标函数和梯度 & 44.4022 & 2571.47 
% & IFM & 8.64566e-07 & 2.292e-05 & 80 & 42.375 \\
% & CCD & 8.99573e-07 & 1.90147e-05 & 700000 & 28.546 \\
% & UCD & 8.48983e-07 & 1.84887e-05 & 990000 & 40.633 \\
% & RCD & 9.45428e-07 & 2.04913e-05 & 980000 & 60.583 \\

% %%%%%%%%%%%%%%%%%%
% 矩阵维数 & 15000 & 12000 
% 最开始的目标函数和梯度 & 50.3221 & 3179.66 
% & IFM & 9.51309e-07 & 2.55215e-05 & 99 & 80.311 \\
% & CCD & 8.7263e-07 & 1.86892e-05 & 1020000 & 45.149 \\
% & UCD & 9.95146e-07 & 2.20648e-05 & 1476000 & 66.916 \\
% & RCD & 9.86487e-07 & 2.20453e-05 & 1440000 & 100.824 \\

%%%%%%%
% 矩阵维数 & 6000 & 600 
% 最开始的目标函数和梯度 & 31.9484 & 449.691 
% & IFM & 28.6202 & 0.00182823 & 40 & 0.587 \\
% & CCD & 28.6202 & 0.000601892 & 8400 & 0.223 \\
% & UCD & 28.6202 & 0.00231571 & 19800 & 0.514 \\
% & RCD & 28.6202 & 0.00134259 & 19800 & 0.55 \\

%%%%%%%%%%%%%%
% 矩阵维数 & 10000 & 2000 
% 最开始的目标函数和梯度 & 40.9439 & 1052.16 
% & IFM & 31.6294 & 0.00295254 & 82 & 7.476 \\
% & CCD & 31.6294 & 0.00166258 & 64000 & 2.348 \\
% & UCD & 31.6294 & 0.00262359 & 136000 & 5.141 \\
% & RCD & 31.6294 & 0.00282412 & 136000 & 5.72 \\

% %%%%%%%%%%%%%
% 矩阵维数 & 15000 & 3000 
% 最开始的目标函数和梯度 & 49.857 & 1570.4 
% & IFM & 38.5489 & 0.00367496 & 86 & 17.019 \\
% & CCD & 38.5489 & 0.00209706 & 102000 & 4.533 \\
% & UCD & 38.5489 & 0.00365486 & 213000 & 9.704 \\
% & RCD & 38.5489 & 0.00388836 & 216000 & 11.137 \\

%%%%%%%%%%%%%%%%
% 矩阵维数 & 20000 & 4000 
% 最开始的目标函数和梯度 & 57.7701 & 2100.37 
% & IFM & 44.9489 & 0.00462223 & 86 & 30.771 \\
% & CCD & 44.9489 & 0.00238819 & 136000 & 7.338 \\
% & UCD & 44.9489 & 0.00485191 & 284000 & 15.694 \\
% & RCD & 44.9489 & 0.00508076 & 280000 & 18.198 \\

%%%%%%%%
% 矩阵维数 & 30000 & 5000 
% 最开始的目标函数和梯度 & 70.5502 & 2888.51 
% & IFM & 57.6083 & 0.00672003 & 67 & 45.666 \\
% & CCD & 57.6083 & 0.00290263 & 130000 & 8.775 \\
% & UCD & 57.6083 & 0.00626242 & 280000 & 19.295 \\
% & RCD & 57.6083 & 0.00609042 & 275000 & 22.521 \\

% %%%%%%%%%%%%%%%%
% 矩阵维数 & 40000 & 6000 
% 最开始的目标函数和梯度 & 81.6042 & 3647 
% & IFM & 68.223 & 0.00833488 & 61 & 65.651 \\
% & CCD & 68.223 & 0.00294496 & 138000 & 11.589 \\
% & UCD & 68.223 & 0.00801527 & 300000 & 24.827 \\
% & RCD & 68.223 & 0.00779252 & 294000 & 28.046 \\
