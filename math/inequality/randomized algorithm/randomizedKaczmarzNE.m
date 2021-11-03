function [x,iter,error_k,iter_k,index_k] = randomizedKaczmarzNE(A, b, x0,maxit,tol,exactx,debug)
%% �����趨
% �������
% A, b, x0 �����ϵ��������ұ��� ��ʼֵ
% maxit,tol,exactx ���������������̶ȣ���ȷ��
% �������
% x iter �������Ľ�, ʵ�ʵ�������
% error ���û�о�ȷ��ʹ洢ÿһ�ε������ݶ�2����
% iter_k ������ݶȵĵ������� Ϊ��ͼ����
% index_k ���������㷨��洢���ѡ�������

%%
[m, n] = size(A);
x = x0;
iter = 0;

%% ����в�
r = b - A * x;
r(r<0) = 0;
norm_Ar = norm(A'*r);
error_k = [norm_Ar];
if ~isempty(exactx)
    e = norm(x-exactx);
    iter_k =[x];
else
    iter_k =[0];
end

index_k=[0];

% ��Ϊ������ֹ������Ҫ����������� Ϊ�˱���ÿ�ε�����ȥ�����ֹ����������ڼ��
iter_test_stop = 1;

alpha = 1;

Arow=sum(A.*A,2);

weight = Arow/sum(Arow);
index=1:m;

for i = 1:maxit
    %randsample to generate weighted random number from given vector
    
    pickedi = randsample(index,1,true,weight);
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


