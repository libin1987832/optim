function [x,iter,error_k,iter_k,index_k] = RFM(A, b, x0, maxit,alpha,maxit_R, tol, exactx,debug)
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
B=A'*A;
colunmnormA=diag(B);
%sAcol = sqrt(colunmnormA);
p=2;
% colunmnormA=power(sqrt(colunmnormA),p);
Ar=A'*r;
for i = 1:maxit
    % ��LSQR�㷨�����������½�����
  %  u = krylovk(A, r, maxit_R);
%     u=wrandomizedGaussSeidel(A, r, Ar, maxit_R,B,colunmnormA,sAcol,p,alpha);
    u=wrandomizedGaussSeidel(A, -r, -Ar, maxit_R,B,colunmnormA,p,alpha);
%    norm(r)
  % [u,~,~]= Rand_Gauss_Seidel_R(A, -r, -Ar,maxit_R,B,colunmnormA,p,alpha);
%     norm(A*u-r)
    % [u,~]= Gass_seidel_D(A, r, maxit_R,colunmnormA,alpha);
    x = x + u;
        r0=r;
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
            iter_k =[iter_k x];
        else
            iter_k =[iter_k i];
        end        
%         error_k = [error_k,e];
            error_k = [error_k,[norm(r0);norm(A*u-r0);norm(r)]];
    end
    
end
end

function x = wrandomizedGaussSeidel(A, r, At_r,maxit,B,Acol,p,alpha)

%%
[m, n] = size(A);


%% ����в�

x = zeros(n,1);

pnormAx_b=(((At_r).^2)./Acol);
prob=(pnormAx_b/sum(pnormAx_b));
cumsumpro=cumsum(prob);
pickedjA=[];
for i = 1:maxit
    rr=rand;
    pickedj=sum(cumsumpro<rr)+1;
   % pickedj=pickedj_a(pickedj_i(i));
  %  pickedj=randsample(index,1,true,weight);
    pickedjA=[pickedjA;pickedj prob' rr];
    col = A(:, pickedj);
    inc = alpha*( col' * r ) / Acol(pickedj);
 %   inc = alpha*( At_r(pickedj) ) / Acol(pickedj);
    x(pickedj) = x(pickedj) - inc;
    r = r - inc*col;
  %  At_r =A'*r;
    At_r = At_r - inc*B(:,pickedj);
    pnormAx_b=(((At_r).^2)./Acol);
%     pnormAx_b=(abs(At_r)./sAcol);
    %pnormAx_b=power((abs(At_r)./sqrt(Acol)),p);
    prob=(pnormAx_b/sum(pnormAx_b));
    cumsumpro=cumsum(prob);
   
end
pickedjA
end

