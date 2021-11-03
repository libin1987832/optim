function [x,iter,error,xA,indexA] = randomizedInexactNE(A, b, x0,maxit,tol,exactx)
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
z=-r;
z(z<0)=0;
z0=z;
r(r<0) = 0;
norm_Ar = norm(A'*r);
error_k = [norm_Ar];
iter_k = [0];
index_k=[0];

% ��Ϊ������ֹ������Ҫ����������� Ϊ�˱���ÿ�ε�����ȥ�����ֹ����������ڼ��
iter_test_stop = 1;
t0=1;


alpha = 1+ min(m/n,n/m);
  Acol=sum(A.*A,1);
  for j = 1:n
      normrow = [normrow,sqrt(Acol(j))];
     index = [index,j];
  end
  weight = normrow/sum(normrow);
if isempty(tol)
  iter = 0;   
  for i = 1:maxit
     pickedj = randsample(index,1,true,weight);
     col = A(:, pickedj);
     rz=(r+z);

     inc = alpha*( col' * rz) / Acol(pickedj);
     x(pickedj) = x(pickedj) + inc;
     r = r - inc*col;
     
     z=-r;
     z(z<0)=0;
     t1 = 0.5 + 0.5 * sqrt(1+4*t0^2);
     z=z0+t0/t1*(z-z0);
     z0=z;
     t0=t1;

    iter =iter+1;
    
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
    %iter = iter+1;
  end

end