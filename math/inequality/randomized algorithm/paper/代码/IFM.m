function [x,iter,error_k,iter_k,index_k] = IFM(A, b, x0, maxit,maxit_LSQR, tol, exactx,debug)
%% �����趨
% �������
% A, b, x0 �����ϵ��������ұ��� ��ʼֵ
% maxit,tol,exactx ���������������̶ȣ���ȷ��
% �������
% x iter �������Ľ�, ʵ�ʵ�������
% error ���û�о�ȷ��ʹ洢ÿһ�ε������ݶ�2����
% iter_k ������ݶȵĵ������� Ϊ��ͼ����
% index_k ���������㷨��洢���ѡ�������

%% �ʼ������
[m,n] = size(A);
x = x0;
iter = 0;

%% �����ʼ�Ĳв�
r = b - A * x0;
r(r<0)=0;
% normAr = norm(A'*r);
norm_r = norm(r);
error_k = [norm_r];
if ~isempty(exactx)
    e = norm(x-exactx);
    iter_k =[x];
else
    iter_k =[0];
end

index_k=[0];

%% �趨��������LSQR�㷨�ĵ�������


Ar=A'*r;
for i = 1:maxit
    % ��LSQR�㷨�����������½�����
    u = krylovk(A, r, maxit_LSQR);
    x = x + u;
    r = b - A * x;
    r( r < 0) = 0;
    iter = iter+1;
    norm_rn = norm(r);
    if abs(norm_rn-norm_r)<tol || norm_rn < tol
%                 fprintf('stop condition:%g,%g',abs(norm_rn-norm_r),norm_rn);
        break;
    end
    norm_r = norm_rn;
 %   if ~isempty(tol) || debug
    if debug
        if ~isempty(exactx)
            e = norm(x-exactx);
            iter_k =[iter_k x];
        else
            iter_k =[iter_k i];
        end        
        error_k = [error_k,e];
    end
    
end
end