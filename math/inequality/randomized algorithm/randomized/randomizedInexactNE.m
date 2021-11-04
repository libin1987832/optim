function [x,iter,error_k,iter_k,index_k] = randomizedInexactNE(A, b, x0,maxit,tol,exactx,debug)
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
z0=-r;
z0(z0<0)=0;
z=z0;
rn=r;
rn(rn<0)=0;

norm_Ar = norm(A'*rn);

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
t0=1;
alpha = 1+ min(m/n,n/m);
%  alpha = 1;
Acol=sum(A.*A,1);
weight = Acol/sum(Acol);
index=1:n;
% pickedj_a =zeros(1,maxit);
% for i = 1:maxit
%  pickedj_a(i) = randsample(index,1,true,weight);
% end
iter_index=1;
for i=1:n-1
    k = ceil(weight(i)*maxit);
    pickedj_a(iter_index:iter_index+k-1)=repmat(i,1,k);
    iter_index=iter_index+k;
end
pickedj_a(iter_index:maxit)=repmat(n,1,maxit-iter_index+1);
pickedj_i=randperm(maxit);
for i = 1:maxit
%     pickedj=pickedj_a(i);
   % pickedj = randsample(index,1,true,weight);
       pickedj=pickedj_a(pickedj_i(i));
    col = A(:, pickedj);
    rz=(r+z);
    
    inc = alpha*( col' * rz) / Acol(pickedj);
    x(pickedj) = x(pickedj) + inc;
    r = r - inc*col;
  %  if mod(iter,1000)==0
   % z=-r;
   % z(z<0)=0;
      z = (-r+abs(r))/2;
    t1 = 0.5 + 0.5 * sqrt(1+4*t0^2);
    z=z0+t0/t1*(z-z0);
    z0=z;
    t0=t1;
 %   end
    
    iter =iter+1;
    
    % ��Ҫ��¼���������е�ֵ ��������
    if mod(iter,iter_test_stop)==0
        if ~isempty(tol) || debug
            rn=r;
            rn(rn<0)=0;
            % �����������
            Ar = A'*rn;
            e = norm(Ar);
            normAr = norm(rn);
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