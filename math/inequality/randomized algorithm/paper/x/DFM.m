function [x,iter,error_k,iter_k,index_k] = DFM(A, b, x0, maxit,alpha,maxit_gs,tol, exactx,debug)
%% �����趨
% Dax ���㷨���� Guass Seidel ���������
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
% error_k = [normAr];
error_k =[];
if ~isempty(exactx)
    e = norm(x-exactx);
    iter_k =[x];
else
    iter_k =[0];
end

index_k=[0];

%% �趨��������LSQR�㷨�ĵ�������
%colunmnormA=sum(A.*A,1);
colunmnormA=[];
for i = 1:n
    colunmnormA = [colunmnormA,norm(A(:,i))^2];
 end
for i = 1:maxit
   u=Gass_seidel_D(A, -r, maxit_gs,colunmnormA,alpha);
    x = x + u;
    r0=r;
    r = b - A * x;
    r( r < 0) = 0;
    norm_rn = norm(r);
    iter = iter+1;
    %terminate condition
    if abs(norm_rn-norm_r)<tol || norm_rn < tol
%                 fprintf('stop condition:%g,%g',abs(norm_rn-norm_r),norm_rn);
        break;
    end
     norm_r = norm_rn;
%     if ~isempty(tol) || debug
    if debug
        if ~isempty(exactx)
            e = norm(x-exactx);
            %iter_k =[iter_k i];
            iter_k =[iter_k x];
        else
            iter_k =[iter_k i];
        end        
          error_k = [error_k,normr];
%         error_k = [error_k,e];
        %    error_k = [error_k,[norm(r0);norm(A*u-r0);norm(r)]];
    end
    
end


end

