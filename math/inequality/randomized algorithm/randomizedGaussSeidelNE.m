function [x,iter,error_k,iter_k,index_k] = randomizedGaussSeidelNE(A, b, x0,maxit,tol,exactx)
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
debug  = 1 ;

%% ����в�
r = b - A * x;
r(r<0) = 0;
norm_Ar = norm(A'*r);
error_k = [norm_Ar];
iter_k = [0];
index_k=[0];
% ��Ϊ������ֹ������Ҫ����������� Ϊ�˱���ÿ�ε�����ȥ�����ֹ����������ڼ��
iter_test_stop = 1;

alpha = 1;

Acol=sum(A.*A,1);

weight = Acol/sum(Acol);


for i = 1:maxit
    
    pickedj = randsample(index,1,true,weight);
    
    col = A(:, pickedj);
    inc = alpha*( col' * r ) / Acol(pickedj);
    x(pickedj) = x(pickedj) + inc;
    r = r - inc*col;
    r( r < 0) = 0;
    iter = iter+1;
    % ��Ҫ��¼���������е�ֵ ��������
    if mod(iter,iter_test_stop)==0
        if ~isempty(tol) || debug
            % ����������� 
            Ar = A'*r;
            e = norm(Ar);
            normAr = norm(r)
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
            end
            error_k = [error_k,e];
            iter_k =[iter_k i];
            index_k = [index_k,pickedj];
        end
    end
end

end