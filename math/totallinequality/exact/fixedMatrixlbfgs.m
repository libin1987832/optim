% 精确方法求解（中间使用MATLAB自带的算法求解带约束的最小二乘问题）
% 将非负最小二乘问题转变成 普通的最小二乘问题?
function [x0, resvec, arvec, face1vec, face2vec, tf]=fixedMatrixlbfgs(A,b,x0,maxit,tol,options)
[m,N]=size(A);
index = 1;
% the residual vector
resvec = zeros(1,maxit + 1);
% the normal gradient 
arvec = zeros(1,maxit + 1);
face1vec = zeros(1,(maxit + 1));
% face x
face2vec = zeros(1,(maxit + 1));
[r, normr, xmin, Ar, KKT, face1, face2] = kktResidual(A, b, x0 , [],1);
resvec(index) = normr;
arvec(index) = KKT;
face1vec(index) = face1;
face2vec(index) = face2;
l  = zeros(N,1);    % lower bound
u  = inf(N,1);      % there is no upper bound

opts    = struct( 'factr', 1e5, 'pgtol', 1e-2, 'm', 50,'verbose',-1);
opts.printEvery     = Inf;
opts.maxIts =400;
t=clock;
AtA     = A'*A; Ab = A'*b;
while KKT > tol && index < maxit
    z = -r;
    z(z<0) = 0;
    bk = b + z;
    Ab = A'*bk;
    %内置使用MATLAB自带的算法
    [x1, ~, info] = lbfgsb(@(x)fungunz(x,A,bk,AtA,Ab), l, u, opts );
  %  [x1,f1,residual,exitflag,output,ff]=lsqlin(A,bk,[],[],[],[],zeros(n,1),Inf*ones(n,1),x0,options);
    x0=x1;
    [r, normr, xmin,Ar, KKT,face1,face2] = kktResidual(A, b, x0 , [], 1);
    %fprintf('index:%d,exit %d,f:%f,res1:%f,res0:%f,ratio:%f!\n',index,exitflag,f1,res1,res0,res1/res0);
    index=index+1;
    resvec(index) = normr;
    arvec(index) = KKT;

    face1vec(index) = face1;
    face2vec(index) = face2;

end
resvec = resvec(resvec>0);
% the normal gradient 
arvec = arvec(arvec>0);
% face b-Ax
face1vec = face1vec(face1vec>0);
% face x
face2vec = face2vec(face2vec>0);

tf = etime(clock,t);
% fprintf('exact run time cost:%f\n',etime(clock,t));
