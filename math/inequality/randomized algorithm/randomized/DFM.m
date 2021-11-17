function [x,iter,error_k,iter_k,index_k] = DFM(A, b, x0, maxit,alpha,maxit_gs,tol, exactx,debug)
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
normAr = norm(A'*r);
error_k = [normAr];
if ~isempty(exactx)
    e = norm(x-exactx);
    iter_k =[x];
else
    iter_k =[0];
end

index_k=[0];

%% �趨��������LSQR�㷨�ĵ�������


% opts.strategy=1;
% opts.p=2;
% 
% opts.lamba = 1+min(n/m,m/n);
B=A'*A;
p=2;
%
colunmnormA=sum(A.*A,1);
Ar=A'*r;
for i = 1:maxit
    % ��LSQR�㷨�����������½�����
  %  u = krylovk(A, r, maxit_R);
   u=Gauss_Seidel(A, -r, maxit_gs,colunmnormA,alpha);
%     norm(r)
 % [u,~]= Gass_seidel_D(A, -r, maxit_gs,colunmnormA,alpha);
%  norm(A*u-r)
% [u,~,~]= Rand_Gauss_Seidel_R(A, r, Ar,maxit_gs,B,colunmnormA',p,alpha);
    x = x + u;
    r = b - A * x;
    r( r < 0) = 0;
    iter = iter+1;
    if ~isempty(tol) || debug
        Ar = A'*r;
        e = norm(Ar);
        normr = norm(r);
        if ~isempty(tol)
            if normr < tol  || e < tol
                break;
            end
        end
    end
    if debug
        if ~isempty(exactx)
            e = norm(x-exactx);
            %iter_k =[iter_k i];
            iter_k =[iter_k x];
        else
            iter_k =[iter_k i];
        end        
        error_k = [error_k,e];
    end
    
end


end
function x= Gauss_Seidel(A, r, maxit,Acol,alpha)

%%
[m, n] = size(A);


%% ����в�
x = zeros(n,1);


for i = 1:maxit
    pickedj=mod(i-1,n)+1;
   % pickedj=pickedj_a(pickedj_i(i));
  %  pickedj=randsample(index,1,true,weight);
    col = A(:, pickedj);
    inc = alpha*( col' * r ) / Acol(pickedj);
  %  inc = alpha*( At_r(pickedj) ) / Acol(pickedj);
    x(pickedj) = x(pickedj) - inc;
    r = r - inc*col;   
end
end