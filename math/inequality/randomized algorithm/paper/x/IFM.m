function [x,iter,error_k,iter_k,index_k] = IFM(A, b, x0, maxit,maxit_LSQR, tol, exactx,debug)
%% 参数设定
% 输入参数
% A, b, x0 问题的系数矩阵和右边项 初始值
% maxit,tol,exactx 最大迭代次数，容忍度，精确解
% 输出参数
% x iter 迭代最后的解, 实际迭代次数
% error 如果没有精确解就存储每一次迭代的梯度2范数
% iter_k 存计算梯度的迭代次数 为作图方便
% index_k 如果是随机算法则存储随机选择的序列

%% 最开始的设置
[m,n] = size(A);
x = x0;
iter = 0;

%% 计算最开始的残差
r = b - A * x0;
r(r<0)=0;
% normAr = norm(A'*r);
norm_r = norm(r);
error_k = [norm_r];
if ~isempty(exactx)
    e = norm(x-exactx);
    iter_k =[x];
else
    iter_k =[0];
end

index_k=[0];

%% 设定子问题中LSQR算法的迭代次数


Ar=A'*r;
for i = 1:maxit
    % 用LSQR算法求解子问题的下降方向
    u = krylovk(A, r, maxit_LSQR);
    x = x + u;
    r = b - A * x;
    r( r < 0) = 0;
    iter = iter+1;
    norm_rn = norm(r);
    if abs(norm_rn-norm_r)<tol || norm_rn < tol
%                 fprintf('stop condition:%g,%g',abs(norm_rn-norm_r),norm_rn);
        break;
    end
    norm_r = norm_rn;
 %   if ~isempty(tol) || debug
    if debug
        if ~isempty(exactx)
            e = norm(x-exactx);
            iter_k =[iter_k x];
        else
            iter_k =[iter_k i];
        end        
        error_k = [error_k,e];
    end
    
end
end