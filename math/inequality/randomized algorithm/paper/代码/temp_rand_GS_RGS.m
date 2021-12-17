clear
clc
debug = 1;
%% 产生问题矩阵
% 随机矩阵
 m = 100;
 n = 10;

 A = 2 * rand(m , n)-1;
% A = [A;-A];
b = 2 * rand(m , 1)-1;
% b = [b;-b];
% b=A*ones(n,1);
x0 = zeros(n , 1);
% save('test2.mat','A','b','x0')
%load('test.mat');
% 二维矩阵
% A = -[1,-1;-1,-1;0,1];b=-[0;-1;0];x0=[-1;0];
% 不一致情况下的正解
% x_exact=[1/2;1/3];

tol=1e-3;
x_exact=[];
%% GuassSeidel
maxit_Rand =500000;
t=clock;
 [x_GS,iter_GS1,error_GS1,xA_GS1,index_GS1] = GuassSeidelNE(A, b, x0,2.0,maxit_Rand,tol,x_exact,debug);
tf_GS=etime(clock,t);
r = b - A * x_GS;
r(r<0) = 0;
r_GS = norm(r);
g_GS = norm(A'*r);
fprintf('& %s & %g & %g & %d & %g \\\\\n', 'CGS(\lambda = 2)', r_GS, g_GS, iter_GS1, tf_GS);


%% simpleGuassSeidel
% maxit_Rand =30000;
t=clock;
[x_GS,iter_GS2,error_GS2,xA_GS2,index_GS2] = simpleGuassSeidelNE(A, b, x0,2,maxit_Rand,tol,x_exact,debug);
tf_GS=etime(clock,t);
r = b - A * x_GS;
r(r<0) = 0;
r_GS = norm(r);
g_GS = norm(A'*r);
fprintf('& %s & %g & %g & %d & %g \\\\\n', 'SGS(\lambda = 2)', r_GS, g_GS, iter_GS2, tf_GS);

%% randGuassSeidel
% maxit_Rand =33000;
t=clock;
[x_GS,iter_GS3,error_GS3,xA_GS3,index_GS3] = randGuassSeidelNE(A, b, x0,2,maxit_Rand,tol,x_exact,debug);
tf_GS=etime(clock,t);
r = b - A * x_GS;
r(r<0) = 0;
r_GS = norm(r);
g_GS = norm(A'*r);
fprintf('& %s & %g & %g & %d & %g \\\\\n', 'RGS(\lambda = 2)', r_GS, g_GS, iter_GS3, tf_GS);


%% 画图
if debug
figure
% h=semilogy(xA_IFM, error_IFM, 'k.');
% h.LineStyle = '--';
 sum = min([iter_GS1,iter_GS2,iter_GS3]) ;
 sum = 20;
 iter_GS1 = sum;
 iter_GS2 = sum;
 iter_GS3 = sum;
display=1:1:iter_GS1;
h=semilogy(xA_GS1(display), error_GS1(display), 'r+');
h.LineStyle = '--';
hold on
display=1:1:iter_GS2;
h=semilogy(xA_GS2(display), error_GS2(display), 'b+');
h.LineStyle = '--';
display=1:1:iter_GS3;
h=semilogy(xA_GS3(display), error_GS3(display), 'g+');
h.LineStyle = '--';
legend('CGS(\lambda = 2)','SGS(\lambda = 2)','RGS(\lambda = 2)');
xlabel('the iterative numbers');
ylabel('the norm of the gradient');
end

error=[error_GS1(display);index_GS1(display);
    error_GS2(display);index_GS2(display);
    error_GS3(display);index_GS3(display)];

error1=[error([1,3,5],2:end) [0;0;0]]-error([1,3,5],1:end);
error=[error;[[0;0;0] error1(:,2:end)]];

