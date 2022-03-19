A=rand(5000,500);
tol=1e-5;
maxit_Rand =2000000;
alpha=0;
error=[];
for i = 1:4
    alpha=alpha+0.5*i;
t=clock;
[x_GS,iter_GS,error_GS,xA_GS,index_GS] = GuassSeidelNE(A, b, x0,alpha,maxit_Rand,tol,x_exact,debug);
tf_GS=etime(clock,t);
r = b - A * x_GS;
r(r<0) = 0;
r_GS = norm(r);
g_GS = norm(A'*r);
fprintf('& %s & %g & %g & %d & %g \\\\\n', 'GuassSeidel', r_GS, g_GS, iter_GS, tf_GS);
error=[error;error_GS];
end
semilogy(error(1,:), 'r') % Matrix A
xlabel('迭代次数') 
ylabel('梯度的误差') 
hold on 
semilogy(error(2,:), 'k') % Matrix M
semilogy(error(3,:), 'g') % Matrix M
semilogy(error(4,:), 'b') % Matrix M
legend('\alpha=0.5', '\alpha=1', '\alpha=1.5', '\alpha=2' )
title('三种指标选择模式对算法的影响')

