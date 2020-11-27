% 精确方法求解（中间使用MATLAB自带的算法求解带约束的最小二乘问题）
% 将非负最小二乘问题转变成 普通的最小二乘问题?
function [x0, resvec, arvec, tf]=fixedMatrix(A,b,x0,maxit,tol,options)
[m,n]=size(A);


index = 1;

% the residual vector
resvec = zeros(1,maxit + 1);
% the normal gradient 
arvec = zeros(1,maxit + 1);


[r, normr, xmin, Ar, KKT] = kktResidual(A, b, x0 , 1);
resvec(index) = normr;
arvec(index) = KKT;
    
t=clock;
while KKT > tol && index < maxit
    z = -r;
    z(z<0) = 0;
    bk = b + z;
    %内置使用MATLAB自带的算法
    [x1,f1,residual,exitflag,output,ff]=lsqlin(A,bk,[],[],[],[],zeros(n,1),Inf*ones(n,1),x0,options);
    %[x1,f1,residual,exitflag]=lsqlin(A,bk,[],[],[],[],[],[],x0,options);
    x0=x1;
    [r, normr, xmin,Ar, KKT] = kktResidual(A, b, x0 ,1);
    %fprintf('index:%d,exit %d,f:%f,res1:%f,res0:%f,ratio:%f!\n',index,exitflag,f1,res1,res0,res1/res0);
    index=index+1;
    resvec(index) = normr;
    arvec(index) = KKT;
end
tf = etime(clock,t);
% fprintf('exact run time cost:%f\n',etime(clock,t));
