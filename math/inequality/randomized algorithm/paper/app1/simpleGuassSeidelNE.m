function [x,iter,error_k,iter_k,index_k] = simpleGuassSeidelNE(A, b, x0,alpha ,maxit,tol,exactx,debug)
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
norm_r = norm(r);
if isempty(exactx)
    error_k = [norm_Ar];
    iter_k =[0];
else
    error_k = [norm(x-exactx)];
    iter_k =[x];
end

index_k=[0];
% ��Ϊ������ֹ������Ҫ����������� Ϊ�˱���ÿ�ε�����ȥ�����ֹ����������ڼ��
iter_test_stop = n;



% Acol=sum(A.*A,1);

Acol = [];
index = [];
% Acol=norm(sum(A.*A,1));
%
% weight = Acol/sum(Acol);
for i = 1:n
    Acol = [Acol,norm(A(:,i))^2];
    index = [index,i];
end

% for i = 1:maxit
%  pickedj_a(i) = randsample(index,1,true,weight);
% end
colrs=0;
for i = 1:maxit
    pickedj = randi(n);
    col = A(:, pickedj);
    colr=col' * r;
    colrs=colrs+colr*colr;
    inc = alpha*( col' * r ) / Acol(pickedj);
    x(pickedj) = x(pickedj) + inc;
    rs = rs - inc*col;
    r = rs;
    %r=(r+abs(r))/2;
    r(r<0)=0;
    iter = iter+1;

    if mod(iter,iter_test_stop)==0
  %     norm_rn = norm(r);
%         if norm_rn > norm_r
%             fprintf('%d',pickedj)
%         end
%         if abs(norm_rn-norm_r)<tol || abs(norm_rn) < tol
%             break;
%         end
%         norm_r = norm_rn;
         if colrs <tol
             break;
         end
         colrs=0;
    end
    % ��Ҫ��¼���������е�ֵ ��������
    if debug
        if mod(iter,iter_test_stop)==0
            if ~isempty(exactx)
                e = norm(x-exactx);
                iter_k =[iter_k x];
            else
                Ar = A'*r;
                e = norm(Ar);
                iter_k =[iter_k i];
            end
            error_k = [error_k,e];
            index_k = [index_k,pickedj];
        end
    end
end
end

