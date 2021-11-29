function [x,iter,error_k,iter_k,index_k] = randomizedGaussSeidelArrNE(A, b, x0,alpha ,p,maxit,tol,exactx,debug)
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
% recurrent Ar
%%
[m, n] = size(A);
x = x0;
iter = 0;

%% 计算残差
r = b - A * x;
rs = r;
r(r<0) = 0;
norm_Ar = norm(A'*r);

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



Acol=sum(A.*A,1)';

At_r=A'*r;
pnormAx_b=power((abs(At_r)./sqrt(Acol)),p);
prob=(pnormAx_b/sum(pnormAx_b));
cumsumpro=cumsum(prob);

N=r>0;
ATNA=A(N,:)'*A(N,:);
for i = 1:maxit
   
    pickedj=sum(cumsumpro<rand)+1;
    
    col = A(:, pickedj);
    inc = alpha*( col' * r ) / Acol(pickedj);
    x(pickedj) = x(pickedj) + inc;
    N= r>0;
    rso=rs;
    rs = rs - inc*col;
    r = rs;
    Nr = r > 0;
    r=(r+abs(r))/2;
    if i >30
    DN=Nr-N;
    sum(abs(DN))
    end
%     if sum(abs(DN))==0
%         At_r = At_r - inc*ATNA(:,pickedj);
%     else
%         DNP = DN>0;
%         if sum(abs(DNP))~=0
%              ADNPT = A(DNP,:)';
%              ATNA = ATNA+ADNPT*ADNPT';
%              At_r = At_r + ADNPT*rso(DNP);
%         end
%         NNP = DN<0;
%         if sum(abs(NNP))~=0
%              ANNPT = A(NNP,:)';
%              ATNA = ATNA-ANNPT*ANNPT';
%              At_r = At_r - ANNPT*rso(NNP);
%         end   
%         At_r = At_r - inc*ATNA(:,pickedj);
%     end
     At_r = A' * r;
   %   end
    pnormAx_b=power((abs(At_r)./sqrt(Acol)),p);
    prob=(pnormAx_b/sum(pnormAx_b));
    cumsumpro=cumsum(prob);
    iter = iter+1;
    % 主要记录迭代过程中的值 用来调试
    if mod(iter,iter_test_stop)==0
        if ~isempty(tol) || debug
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
            index_k = [index_k,pickedj];
        end
    end
end

end