function [x,iter,error_k,iter_k,index_k] = RFM(A, b, x0, maxit,alpha,maxit_R, tol, exactx,debug)
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
normAr = norm(A'*r);
error_k = [normAr];
if ~isempty(exactx)
    e = norm(x-exactx);
    iter_k =[x];
else
    iter_k =[0];
end

index_k=[0];

%% 设定子问题中LSQR算法的迭代次数
B=A'*A;
colunmnormA=diag(B);
p=2;
% colunmnormA=power(sqrt(colunmnormA),p);
Ar=A'*r;
for i = 1:maxit
    % 用LSQR算法求解子问题的下降方向
  %  u = krylovk(A, r, maxit_R);
    u=wrandomizedGaussSeidel(A, r, Ar, maxit_R,B,colunmnormA,p,alpha);
    x = x + u;
    r = b - A * x;
    r( r < 0) = 0;
    iter = iter+1;
    if ~isempty(tol) || debug
        Ar = A'*r;
        e = norm(Ar);
        normr = norm(r);
        if ~isempty(tol)
            if normr < tol  || e < tol
                break;
            end
        end
    end
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

function x = wrandomizedGaussSeidel(A, r, At_r,maxit,B,Acol,p,alpha)

%%
[m, n] = size(A);


%% 计算残差

x = zeros(n,1);

%pnormAx_b=power((abs(At_r)./sqrt(Acol)),p);

%pnormAx_b=((At_r.^p)./Acol);
prob=(pnormAx_b/sum(pnormAx_b));
 pnormAx_b=(abs(At_r)./sAcol);
cumsumpro=cumsum(prob);

for i = 1:maxit
    pickedj=sum(cumsumpro<rand)+1;
   % pickedj=pickedj_a(pickedj_i(i));
  %  pickedj=randsample(index,1,true,weight);
    
    col = A(:, pickedj);
    inc = alpha*( col' * r ) / Acol(pickedj);
 %   inc = alpha*( At_r(pickedj) ) / Acol(pickedj);
    x(pickedj) = x(pickedj) + inc;
    r = r - inc*col;
    
    At_r = At_r - inc*B(:,pickedj);
    pnormAx_b=(abs(At_r)./sAcol);
    %pnormAx_b=power((abs(At_r)./sqrt(Acol)),p);
    prob=(pnormAx_b/sum(pnormAx_b));
    cumsumpro=cumsum(prob);
   
end

end

