function [x,iter,error_k,iter_k,index_k] = swrandKacz(A, b, x0,maxit,tol,exactx,debug)
% || (b-Ax)_+ || search weight
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

alpha = 1;

AAT=A*A';
Arow=sum(A.*A,2);

weight = Arow/sum(Arow);
index=1:m;

Ip = A>0;
Im = A<0;
Ie = abs(A)<1e-15;
I  = eye(n);

for i = 1:maxit
    col = A(:, pickedj);
    Icolp = Ip(:,pickedj);
    Icolm = Im(:,pickedj);
    palpha=rs(Icolp)./col(Icolp);
    maxrcol = max(palpha);
    malpha = rs(Icolm)./col(Icolm);
    minrcol = min(malpha);

    if sum(Icolm)==0
        inc = maxrcol;
    elseif sum(Icolp)==0
        inc = minrcol;
    elseif maxrcol < minrcol
        inc = maxrcol;
    else
       active1 = plaha(minrcol < plaha);
       active2 = malpha(malpha < maxrcol);
       active = sort([active1,active2],'ascend');
       AN_r+r(i)*A(i,:)' - active*A(:,i)'*col
      %  inc = spiecewise(A,b,s*I(:,pickedj),x);
      inc = bisect2(minrcol,maxrcol,A,b,x,I(:,pickedj),1e-10);
    end
    
    
    row = A(pickedi, :);
    r_pickedi=b(pickedi) - (row * x);
    if r_pickedi>0
        x = x + alpha * ( r_pickedi ) / (Arow(pickedi)) * row';
    end
    iter = iter+1;
    
    % ��Ҫ��¼���������е�ֵ ��������
    if mod(iter,iter_test_stop)==0
        if ~isempty(tol) || debug
            r = b - A * x;
            r(r<0) = 0;
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
            index_k = [index_k,pickedi];
        end
    end
    
end


