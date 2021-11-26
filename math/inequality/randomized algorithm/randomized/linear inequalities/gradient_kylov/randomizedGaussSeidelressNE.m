function [x,iter,error_k,iter_k,index_k] = randomizedGaussSeidelressNE(A, b, x0,alpha ,maxit,tol,exactx,debug)
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

p=2;
r_orig = r;
x_orig = x;
rs_orig = rs;
normar = zeros(n,1);
for pickedj_f = 1:n
    col = A(:, pickedj_f);
    inc = alpha*( col' * r ) / Acol(pickedj_f);
    x(pickedj_f) = x(pickedj_f) + inc;
    rs = rs - inc*col;
    r = rs;
    r=(r+abs(r))/2;
    normar(pickedj_f) = norm(r);
    r= r_orig;
    x = x_orig;
    rs = rs_orig;
end
normar = norm(r_orig)-normar;
prob=(normar/sum(normar));
cumsumpro=cumsum(prob);
iter_index = 1;
for i=1:n-1
    k = ceil(maxit/n);
    pickedj_a(iter_index:iter_index+k-1)=repmat(i,1,k);
    iter_index=iter_index+k;
end
pickedj_a(iter_index:maxit)=repmat(n,1,maxit-iter_index+1);
%pickedj_i=randperm(maxit);
k=ceil(maxit/n);
for i = 1:maxit
 %   pickedj=sum(cumsumpro<rand)+1;
 mod(i,n)*k+ceil(i/n)
   pickedj=pickedj_a((mod(i,n))*k+mod(i,n));
    col = A(:, pickedj);
    colr = col' * r;
    inc = alpha*( colr ) / Acol(pickedj);
    x(pickedj) = x(pickedj) + inc;
    rs = rs - inc*col;
    normor=norm(r);
    r = rs;
    r=(r+abs(r))/2;
    
    normar(pickedj) = normor-norm(r);

    prob=(normar/sum(normar));
    cumsumpro=cumsum(prob);
    %   end
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