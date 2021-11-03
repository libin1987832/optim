function [x,iter,error,xA,indexA] = randomizedInexactNE(A, b, x0,maxit,tol,exactx)
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
debug  = 1 ;

%% 计算残差
r = b - A * x;
z=-r;
z(z<0)=0;
z0=z;
r(r<0) = 0;
norm_Ar = norm(A'*r);
error_k = [norm_Ar];
iter_k = [0];
index_k=[0];

% 因为测试终止条件需要矩阵乘以向量 为了避免每次迭代都去检测终止条件因此周期检测
iter_test_stop = 1;
t0=1;


alpha = 1+ min(m/n,n/m);
  Acol=sum(A.*A,1);
  for j = 1:n
      normrow = [normrow,sqrt(Acol(j))];
     index = [index,j];
  end
  weight = normrow/sum(normrow);
if isempty(tol)
  iter = 0;   
  for i = 1:maxit
     pickedj = randsample(index,1,true,weight);
     col = A(:, pickedj);
     rz=(r+z);

     inc = alpha*( col' * rz) / Acol(pickedj);
     x(pickedj) = x(pickedj) + inc;
     r = r - inc*col;
     
     z=-r;
     z(z<0)=0;
     t1 = 0.5 + 0.5 * sqrt(1+4*t0^2);
     z=z0+t0/t1*(z-z0);
     z0=z;
     t0=t1;

    iter =iter+1;
    
    % 主要记录迭代过程中的值 用来调试
    if mod(iter,iter_test_stop)==0
        if ~isempty(tol) || debug
            % 矩阵乘以向量 
            Ar = A'*r;
            e = norm(Ar);
            normAr = norm(r)
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
            end
            error_k = [error_k,e];
            iter_k =[iter_k i];
            index_k = [index_k,pickedj];
        end
    end
    %iter = iter+1;
  end

end