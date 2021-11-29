function [x,iter,error_k,iter_k,index_k] = optimGaussSeidel(A, b, x0,alpha ,maxit,tol,exactx,debug)
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
%index_k=[zeros(n+1,1)];
% ��Ϊ������ֹ������Ҫ����������� Ϊ�˱���ÿ�ε�����ȥ�����ֹ����������ڼ��
iter_test_stop = 1;



Acol=sum(A.*A,1);

weight = Acol/sum(Acol);
index=1:n;
pickedj_a =zeros(1,maxit);
% for i = 1:maxit
%  pickedj_a(i) = randsample(index,1,true,weight);
% end
iter_index=1;
normr_opt=norm(r);
for i = 1:maxit
    pickedj = 1;
    r_orig = r;
    x_orig = x;
    rs_orig = rs;
  %  normar = zeros(n,1);
    for pickedj_f = 1:n
        col = A(:, pickedj_f);
        inc = alpha*( col' * r ) / Acol(pickedj_f);
        x(pickedj_f) = x(pickedj_f) + inc;
        rs = rs - inc*col;
        r = rs;
        r=(r+abs(r))/2;
        normr = norm(r);
       % normar(pickedj_f) = normr;
        if normr<normr_opt
            pickedj = pickedj_f;
            r_opt = r;
            x_opt = x;
            rs_opt = rs;
            normr_opt = normr;
        end
        r= r_orig;
        x = x_orig;
        rs = rs_orig;
    end
    r = r_opt;
    x = x_opt;
    rs = rs_opt;
   %   end
    iter = iter+1;
    % ��Ҫ��¼���������е�ֵ ��������
    if mod(iter,iter_test_stop)==0
        if ~isempty(tol) || debug
            % ����������� 
            Ar = A'*r;
            e = norm(Ar);
            normAr = norm(r);
            e = normAr;
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
         %   index_k = [index_k,[pickedj;normar]];
        end
    end
end

end