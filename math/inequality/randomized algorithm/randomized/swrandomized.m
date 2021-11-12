function [x,iter,error_k,iter_k,index_k] = swrandomized(A, b, x0,p,maxit,tol,exactx,debug)
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

% weight = Acol/sum(Acol);
% index=1:n;
% pickedj_a =zeros(1,maxit);
% for i = 1:maxit
%  pickedj_a(i) = randsample(index,1,true,weight);
% end
% iter_index=1;
At_r=A'*r;
pnormAx_b=power((abs(At_r)./sqrt(Acol)),p);
prob=(pnormAx_b/sum(pnormAx_b));
cumsumpro=cumsum(prob);

Ip = A>0;
Im = A<0;
Ie = abs(A)<1e-15;
I  = eye(n);
for i = 1:maxit
    %pickedj=pickedj_a(pickedj_i(i));
    pickedj=sum(cumsumpro<rand)+1;
  %  pickedj=randsample(index,1,true,weight);
    
    col = A(:, pickedj);
    Icolp = Ip(:,pickedj);
    Icolm = Im(:,pickedj);
    maxrcol = max(rs(Icolp)./col(Icolp));
    minrcol = min(rs(Icolm)./col(Icolm));

    if sum(Icolm)==0
        inc = maxrcol;
    elseif sum(Icolp)==0
        inc = minrcol;
    elseif maxrcol < minrcol
        inc = maxrcol;
    else
      %  inc = spiecewise(A,b,s*I(:,pickedj),x);
      inc = bisect2(minrcol,maxrcol,A,b,x,I(:,pickedj),1e-10);
    end

   % inc = alpha*( col' * r ) / Acol(pickedj);
   


    x(pickedj) = x(pickedj) + inc;
    rs = rs - inc*col;
    

     r = rs;
   %    if mod(iter,100)==0
       % r( r < 0) = 0;
       r=(r+abs(r))/2;
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
function xc=bisect2(a,b,A,Ab,x0,p,tol)
if sign(f(A,Ab,x0,p,a))*sign(f(A,Ab,x0,p,b)) >= 0
    error('f(a)f(b)<0 not satisfied!') %ceases execution
end
fa=f(A,Ab,x0,p,a);
fb=f(A,Ab,x0,p,b);
while (b-a)/2>tol
    c=(a+b)/2;
    fc=f(A,Ab,x0,p,c);
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
function fx = f(A,b,x0,p,x)
    x=x0+x*p;
    r= b-A*x;
    r(r<0)=0;
    fx = p'*A'*r;
end
end