function [x,iter,error_k,iter_k,index_k] = DFM(A, b, x0, maxit,alpha,maxit_gs,tol, exactx,debug)
%% 参数设定
% Dax 的算法采用 Guass Seidel 求解子问题
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
% error_k = [normAr];
error_k =[];
if ~isempty(exactx)
    e = norm(x-exactx);
    iter_k =[x];
else
    iter_k =[0];
end

index_k=[0];

%% 设定子问题中LSQR算法的迭代次数
%colunmnormA=sum(A.*A,1);
colunmnormA=[];
for i = 1:n
    colunmnormA = [colunmnormA,norm(A(:,i))^2];
 end
for i = 1:maxit
   u=Gass_seidel_D(A, -r, maxit_gs,colunmnormA,alpha);
    x = x + u;
    r0=r;
    r = b - A * x;
    r( r < 0) = 0;
    norm_rn = norm(r);
    iter = iter+1;
    %terminate condition
    if abs(norm_rn-norm_r)<tol || norm_rn < tol
%                 fprintf('stop condition:%g,%g',abs(norm_rn-norm_r),norm_rn);
        break;
    end
     norm_r = norm_rn;
%     if ~isempty(tol) || debug
    if debug
        if ~isempty(exactx)
            e = norm(x-exactx);
            %iter_k =[iter_k i];
            iter_k =[iter_k x];
        else
            iter_k =[iter_k i];
        end        
          error_k = [error_k,normr];
%         error_k = [error_k,e];
        %    error_k = [error_k,[norm(r0);norm(A*u-r0);norm(r)]];
    end
    
end


end

