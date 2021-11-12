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
    s=sign(col'*r);
    if minrcol <0 
        minrcol = 0;
    end
    if maxrcol <0
        maxrcol=0;
    end
    if sum(Icolm)==0
        inc = maxrcol;
    end
    if sum(Icolp)==0
        inc = minrcol;
    end
    if maxrcol < minrcol
        inc = maxrcol;
    else
        inc = spiecewise(A,b,s*I(:,pickedj),x);
    end

   % inc = alpha*( col' * r ) / Acol(pickedj);
   


    x(pickedj) = x(pickedj) + s*inc;
    rs = rs - s*inc*col;
    

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

end