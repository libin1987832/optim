function [x,iter,error_k,iter_k,index_k] = randomizedKaczmarzNE(A, b, x0,maxit,tol,exactx,debug)
%% 参数设定
% 输入参数
% A, b, x0 问题的系数矩阵和右边项 初始值
% maxit,tol,exactx 最大迭代次数，容忍度，精确解
% 输出参数
% x iter 迭代最后的解, 实际迭代次数
% error 如果没有精确解就存储每一次迭代的梯度2范数
% iter_k 存计算梯度的迭代次数 为作图方便
% index_k 如果是随机算法则存储随机选择的序列

%%
[m, n] = size(A);
x = x0;
iter = 0;

%% 计算残差
r = b - A * x;
r(r<0) = 0;
norm_Ar = norm(A'*r);
error_k = [norm_Ar];
if ~isempty(exactx)
    e = norm(x-exactx);
    iter_k =[x];
else
    iter_k =[0];
end

index_k=[0];

% 因为测试终止条件需要矩阵乘以向量 为了避免每次迭代都去检测终止条件因此周期检测
iter_test_stop = 1;

alpha = 1;

Arow=sum(A.*A,2);

weight = Arow/sum(Arow);
index=1:m;

for i = 1:maxit
    %randsample to generate weighted random number from given vector
    
    pickedi = randsample(index,1,true,weight);
    row = A(pickedi, :);
    r_pickedi=b(pickedi) - (row * x);
    if r_pickedi>0
        x = x + alpha * ( r_pickedi ) / (Arow(pickedi)) * row';
    end
    iter = iter+1;
    
    % 主要记录迭代过程中的值 用来调试
    if mod(iter,iter_test_stop)==0
        if ~isempty(tol) || debug
            r = b - A * x;
            r(r<0) = 0;
            % 矩阵乘以向量 
            Ar = A'*r;
            e = norm(Ar);
            normAr = norm(r);
            % 如果有容忍度 即使没有到达最大迭代次数也终止
            if ~isempty(tol)
                if normAr < tol  || e < tol
                    break;
                end
            end
        end
        % 用来记录迭代信息 可能会影响效率
        if debug
            if ~isempty(exactx)
                e = norm(x-exactx);
                iter_k =[iter_k x];
            else
                iter_k =[iter_k i];
            end
            error_k = [error_k,e];
            index_k = [index_k,pickedi];
        end
    end
    
end


