function [x,iter,error_k,iter_k,index_k] = randomizedInexactNE(A, b, x0,maxit,tol,exactx,debug)
%% 参数设定
% 输入参数
% A, b, x0 问题的系数矩阵和右边项 初始值
% maxit,tol,exactx 最大迭代次数，容忍度，精确解
% 输出参数
% x iter 迭代最后的解, 实际迭代次数
% error 如果没有精确解就存储每一次迭代的梯度2范数
% iter_k 存计算梯度的迭代次数 为作图方便
% index_k 如果是随机算法则存储随机选择的序列
% exactx 精确值为空则将iter_k 放每次迭代的x值
% debug 调试程序
%%
[m, n] = size(A);
x = x0;
iter = 0;


%% 计算残差
r = b - A * x;
z0=-r;
z0(z0<0)=0;
z=z0;
rn=r;
rn(rn<0)=0;

norm_Ar = norm(A'*rn);

if isempty(exactx)
    error_k = [norm_Ar];
    iter_k =[0];
else
    error_k = [norm(x-exactx)];
    iter_k =[x];
end

index_k=[0];

% 因为测试终止条件需要矩阵乘以向量 为了避免每次迭代都去检测终止条件因此周期检测
iter_test_stop = 1;
t0=1;
alpha = 1+ min(m/n,n/m);
%  alpha = 1;
Acol=sum(A.*A,1);
weight = Acol/sum(Acol);
index=1:n;
% pickedj_a =zeros(1,maxit);
% for i = 1:maxit
%  pickedj_a(i) = randsample(index,1,true,weight);
% end
iter_index=1;
for i=1:n-1
    k = ceil(weight(i)*maxit);
    pickedj_a(iter_index:iter_index+k-1)=repmat(i,1,k);
    iter_index=iter_index+k;
end
pickedj_a(iter_index:maxit)=repmat(n,1,maxit-iter_index+1);
pickedj_i=randperm(maxit);
for i = 1:maxit
%     pickedj=pickedj_a(i);
   % pickedj = randsample(index,1,true,weight);
       pickedj=pickedj_a(pickedj_i(i));
    col = A(:, pickedj);
    rz=(r+z);
    
    inc = alpha*( col' * rz) / Acol(pickedj);
    x(pickedj) = x(pickedj) + inc;
    r = r - inc*col;
  %  if mod(iter,1000)==0
   % z=-r;
   % z(z<0)=0;
      z = (-r+abs(r))/2;
    t1 = 0.5 + 0.5 * sqrt(1+4*t0^2);
    z=z0+t0/t1*(z-z0);
    z0=z;
    t0=t1;
 %   end
    
    iter =iter+1;
    
    % 主要记录迭代过程中的值 用来调试
    if mod(iter,iter_test_stop)==0
        if ~isempty(tol) || debug
            rn=r;
            rn(rn<0)=0;
            % 矩阵乘以向量
            Ar = A'*rn;
            e = norm(Ar);
            normAr = norm(rn);
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
            index_k = [index_k,pickedj];
        end
    end
end
