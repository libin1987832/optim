function [x,iter,error_k,iter_k,index_k] = randomizedGaussSeidelArrNE(A, b, x0,alpha ,p,maxit,tol,exactx,debug)
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
% recurrent Ar
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

At_r=A'*r;
pnormAx_b=power((abs(At_r)./sqrt(Acol)),p);
prob=(pnormAx_b/sum(pnormAx_b));
cumsumpro=cumsum(prob);

N=r>0;
ATNA=A(N,:)'*A(N,:);
for i = 1:maxit
   
    pickedj=sum(cumsumpro<rand)+1;
    
    col = A(:, pickedj);
    inc = alpha*( col' * r ) / Acol(pickedj);
    x(pickedj) = x(pickedj) + inc;
    N= r>0;
    rso=rs;
    rs = rs - inc*col;
    r = rs;
    Nr = r > 0;
    r=(r+abs(r))/2;
    if i >30
    DN=Nr-N;
    sum(abs(DN))
    end
%     if sum(abs(DN))==0
%         At_r = At_r - inc*ATNA(:,pickedj);
%     else
%         DNP = DN>0;
%         if sum(abs(DNP))~=0
%              ADNPT = A(DNP,:)';
%              ATNA = ATNA+ADNPT*ADNPT';
%              At_r = At_r + ADNPT*rso(DNP);
%         end
%         NNP = DN<0;
%         if sum(abs(NNP))~=0
%              ANNPT = A(NNP,:)';
%              ATNA = ATNA-ANNPT*ANNPT';
%              At_r = At_r - ANNPT*rso(NNP);
%         end   
%         At_r = At_r - inc*ATNA(:,pickedj);
%     end
     At_r = A' * r;
   %   end
    pnormAx_b=power((abs(At_r)./sqrt(Acol)),p);
    prob=(pnormAx_b/sum(pnormAx_b));
    cumsumpro=cumsum(prob);
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