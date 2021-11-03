clear
clc
%% 产生问题矩阵


% 二维矩阵
 A = -[1,-1;-1,-1;0,1];b=-[0;-1;0];x0=[-1;0];
% 不一致情况下的正解
% x_exact=[1/2;1/3];
% 一致情况下的正解
% x_exact=[0;0];

%% 基于IFM的算法找到一个解
[~, ~, norm_r, norm_g] = residual(A,b,x0);
fprintf('%s & %g & %g \n','最开始的目标函数和梯度', norm_r, norm_g);
[x_exact, ~, ~, ~, ~] = IFM(A, b, x0,1000,1e-15,[]);
[~, ~, norm_r, norm_g] = residual(A,b,x_exact);
fprintf('%s & %g & %g \n','IFM解的目标函数值和梯度  ', norm_r, norm_g);

%% 参数的设定
maxit_IFM = 100;
maxit_Rand = 100;
tol=1e-5;

%% IFM算法求解问题
t=clock;
[x_IFM,iter_IFM,error_IFM,xA_IFM,index_IFM] = IFM(A, b, x0, maxit_IFM, tol, x_exact);
tf_IFM=etime(clock,t);
[~, ~, r_IFM, g_IFM] = residual(A,b,x_IFM);
fprintf('& %s & %s & %s & %s & %s \\\\\n', 'alg', 'norm(r_+)', 'norm(Ar_+)', 'iteration', 'time');
fprintf('& %s & %g & %g & %d & %g \\\\\n','IFM', r_IFM, g_IFM,iter_IFM,tf_IFM);

%% Kaczmarz
t=clock;
[x_Kac,iter_Kac,error_Kac,xA_Kac,index_Kac] = randomizedKaczmarzNE(A, b, x0, maxit_Rand,tol,x_exact);
tf_Kac=etime(clock,t);
[~, ~, r_Kac, g_Kac] = residual(A, b, x_Kac);
fprintf('& %s & %g & %g & %d & %g \\\\\n', 'Kaczmarz', r_Kac, g_Kac, iter_Kac, tf_Kac);

%% GaussSeidel
t=clock;
[x_GS,iter_GS,error_GS,xA_GS,index_GS] = randomizedGaussSeidelNE(A, b, x0, maxit_Rand, tol,x_exact);
tf_GS=etime(clock,t);
[~, ~, r_GS, g_GS] = residual(A, b, x_GS);
fprintf('& %s & %g & %g & %d & %g \\\\\n', 'Gauss', r_GS, g_GS, iter_GS, tf_GS);

%% InexactGaussSeidel
t=clock;
[x_In,iter_In,error_In,xA_In,index_In] = randomizedInexactNE(A, b, x0,maxit_Rand,tol,x_exact);
tf_In=etime(clock,t);
[~, ~, r_In, g_In] = residual(A, b, x_In);
fprintf('& %s & %g & %g & %d & %g \\\\\n','Inexact', r_In, g_In,iter_In, tf_In);


x = linspace(-1,1.5);
y = linspace(-0.5,1);
[X,Y] = meshgrid(x,y);
XZ = repmat(X,1,1,3);
YZ = repmat(Y,1,1,3);
ba1 = reshape(A(:,1),1,1,3);
ba2 = reshape(A(:,2),1,1,3);
br = reshape(b,1,1,3);
z = bsxfun(@times,ba1,XZ )+bsxfun(@times,ba2,YZ );
z = bsxfun(@minus,br,z);
z(z<0)=0;
z=0.5*z.^2;
Z=squeeze(sum(z,3));
figure
contour(X,Y,Z,40)
hold on
plot(xA_Kac(1,:),xA_Kac(2,:),'b+')
% plot(xA_GS(1,:),xA_GS(2,:),'ro')
%plot(xA_In(1,:),xA_In(2,:),'g*')
line([-0.5,1],[-0.5,1]);
line([0,1.5],[1,-0.5]);
line([-1,1.5],[0,0]);

