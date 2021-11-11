function [x,iter,error_k,iter_k,index_k] = wrandomizedGaussSeidel(A, b, x0,p,maxit,tol,exactx,debug)
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
r = A * x - b;
norm_Ar = norm(r);
r0=norm_Ar;
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

B = A'*A;
Acol = diag(B);
% Acol=sum(A.*A,1);

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

for i = 1:maxit
    pickedj=sum(cumsumpro<rand)+1;
   % pickedj=pickedj_a(pickedj_i(i));
  %  pickedj=randsample(index,1,true,weight);
    
  col = A(:, pickedj);
%    inc = alpha*( col' * r ) / Acol(pickedj
    inc = alpha*( At_r(pickedj) ) / Acol(pickedj);
    x(pickedj) = x(pickedj) - inc;
    r = r - inc*col;
    
    At_r = At_r - inc*B(:,pickedj);
    pnormAx_b=power((abs(At_r)./sqrt(Acol)),p);
    prob=(pnormAx_b/sum(pnormAx_b));
    cumsumpro=cumsum(prob);
    
    iter = iter+1;
    % ��Ҫ��¼���������е�ֵ ��������
    if mod(iter,iter_test_stop)==0
        if ~isempty(tol) || debug
            % ����������� 
            e = norm(r)/r0;
            % ��������̶� ��ʹû�е�������������Ҳ��ֹ
            if ~isempty(tol)
                if e < tol
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