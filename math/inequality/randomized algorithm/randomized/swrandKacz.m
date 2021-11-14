function [x,iter,error_k,iter_k,index_k] = swrandKacz(A, b, x0,p,maxit,tol,exactx,debug)
% || (b-Ax)_+ || search weight
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
rs = r;
I_old = r<0;
r(I_old ) = 0;
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


AAT=A*A';

index=1:m;


% Ie = abs(A)<1e-15;
% I  = eye(n);
   pnormAx_b=power((abs(r)./sqrt(r)),p);
    prob=(pnormAx_b/sum(pnormAx_b));
    cumsumpro=cumsum(prob);
for i = 1:maxit
    pickedi=sum(cumsumpro<rand)+1;
    AATrow = AAT(:, pickedi);
    Icolp = AATrow>0;
    Icolm =  AATrow<0;
    palpha=rs(Icolp)./AATrow(Icolp);
    maxrcol = max(palpha);
    malpha = rs(Icolm)./AATrow(Icolm);
    minrcol = min(malpha);

    if sum(Icolm)==0
        inc = maxrcol;
    elseif sum(Icolp)==0
        inc = minrcol;
    elseif maxrcol < minrcol
        inc = maxrcol;
    else
           alphSort=sort([minrcol;palpha(minrcol < palpha);malpha(malpha < maxrcol);maxrcol]);
       for j = 1:size(alphSort,1)
           rj=b-alphSort(j)*AATrow;
           df=AATrow'*rj;
           if df >0
               inc = 0.5*(alphSort(j)+alphSort(j-1));
                break;
           end
       end
      %  inc = spiecewise(A,b,s*I(:,pickedj),x);
%       inc = bisect2(minrcol,maxrcol,rs,AATrow,1e-1);
    end
    
    
    row = A(pickedi, :);
    x = x + inc * row';
    rs = rs - inc * AATrow;
    r=rs;
    I=r<0;
    r(I)=0;
    %I=rs<0;
    pnormAx_b=power((abs(r)./sqrt(r)),p);
    prob=(pnormAx_b/sum(pnormAx_b));
    cumsumpro=cumsum(prob);
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
function xc=bisect2(a,b,r,row,tol)
fa=f(r,row,a);
fb=f(r,row,b);    
if sign(fa)*sign(fb) >= 0
    error('f(a)f(b)<0 not satisfied!') %ceases execution
end
while (b-a)/2>tol
    c=(a+b)/2;
    fc=f(r,row,c);
    if fc == 0              %c is a solution, done
        break
    end
    if sign(fc)*sign(fa)<0  %a and c make the new interval
        b=c;
        fb=fc;
    else%c and b make the new interval
        a=c;
        fa=fc;
    end
end
xc=(a+b)/2;
end
function fx = f(r,row,x)
    r= r-x*row;
    r(r<0)=0;
    fx = row'*r;
end

end