function [x,iter,error_k,iter_k,index_k] = swrandomized(A, b, x0,p,maxit,tol,exactx,debug)
%% �����趨
% �������
% A, b, x0 �����ϵ��������ұ��� ��ʼֵ
% maxit,tol,exactx ���������������̶ȣ���ȷ��
% �������
% x iter �������Ľ�, ʵ�ʵ�������
% error ���û�о�ȷ��ʹ洢ÿһ�ε������ݶ�2����
% iter_k ������ݶȵĵ������� Ϊ��ͼ����
% index_k ���������㷨��洢���ѡ�������
% exactx ��ȷֵΪ����iter_k ��ÿ�ε�����xֵ
% debug ���Գ���
%%
[m, n] = size(A);
x = x0;
iter = 0;

%% ����в�
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
% ��Ϊ������ֹ������Ҫ����������� Ϊ�˱���ÿ�ε�����ȥ�����ֹ����������ڼ��
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
    % ��Ҫ��¼���������е�ֵ ��������
    if mod(iter,iter_test_stop)==0
        if ~isempty(tol) || debug
            % ����������� 
            Ar = A'*r;
            e = norm(Ar);
            normAr = norm(r);
            % ��������̶� ��ʹû�е�������������Ҳ��ֹ
            if ~isempty(tol)
                if normAr < tol  || e < tol
                    break;
                end
            end
        end
        % ������¼������Ϣ ���ܻ�Ӱ��Ч��
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