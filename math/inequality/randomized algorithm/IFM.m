function [x,iter,error_k,iter_k,index_k] = IFM(A, b, x0, maxit, tol, exactx)
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
debug = 0;
%% �����ʼ�Ĳв�
r = b - A * x0;
r(r<0)=0;
normAr = norm(A'*r);
error_k = [normAr];
iter_k = [0];
index_k=[0];

%% �趨��������LSQR�㷨�ĵ�������
maxit_LSQR=3;


for i = 1:maxit
    % ��LSQR�㷨�����������½�����
    u = krylovk(A, r, maxit_LSQR);
    x = x + u;
    r = b - A * x;
    r( r < 0) = 0;
    iter = iter+1;
    if ~isempty(tol) || debug
        Ar = A'*r;
        e = norm(Ar);
        normAr = norm(r)
        if ~isempty(tol) 
            if normAr < tol  || e < tol
              break;
            end
        end
    end
    if debug
        if ~isempty(exactx)
            e = norm(x-exactx);
        end
        error_k = [error_k,e];
        iter_k =[iter_k i];
    end
    
end
% ��LSQR�㷨���������� argmin F��u�� = || Au + Ax_k-b-z_k || = || Au -y_k || 

    function xk=krylovk(A,y,k)
        u1=0;
        
        beta1=norm(y);q1=y/beta1;v1=A'*q1;alph1=norm(v1);v1=v1/alph1;
        
        ro_1=alph1;
        thgma_1=beta1;
        g1=v1;
        
        for i=1:k
            q2=A*v1-alph1*q1;beta2=norm(q2);q2=q2/beta2;
            
            ro1=norm([ro_1,beta2]);c1=ro_1/ro1;s1=beta2/ro1;
            
            v2=A'*q2-beta2*v1;alph2=norm(v2);v2=v2/alph2;
            
            theta2=s1*alph2;ro_2=c1*alph2;thgma1=c1*thgma_1;thgma_2=-1*s1*thgma_1;
            
            u2=u1+thgma1*g1./ro1;  g2=v2-theta2*g1./ro1;
            
            u1=u2;
            q1=q2;v1=v2;alph1=alph2;
            
            ro_1=ro_2;
            thgma_1=thgma_2;
            g1=g2;
        end
        xk=u1;
    end

end