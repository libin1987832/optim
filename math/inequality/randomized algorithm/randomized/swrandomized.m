function [x,iter,error_k,iter_k,index_k] = swrandomized(A, b, x0,p,alpha,maxit,tol,exactx,debug)
% search weight Guass 
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
I_old= r>0;
r=0.5*(r+abs(r));
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


ATA=A'*A;
Acol=diag(ATA);
%Acol=sum(A.*A,1)';

% weight = Acol/sum(Acol);
% index=1:n;
% pickedj_a =zeros(1,maxit);
% for i = 1:maxit
%  pickedj_a(i) = randsample(index,1,true,weight);
% end
% iter_index=1;
At_r=A'*r;
pnormAx_b=power((abs(At_r)./sqrt(Acol)),p);
prob=(pnormAx_b/sum(pnormAx_b));
cumsumpro=cumsum(prob);

Ip = A>0;
Im = A<0;
Ie = abs(A)<1e-15;
I  = eye(n);
for i = 1:maxit
    %pickedj=pickedj_a(pickedj_i(i));
    pickedj=sum(cumsumpro<rand)+1;
  %  pickedj=randsample(index,1,true,weight);
    
    col = A(:, pickedj);
%     Icolp = Ip(:,pickedj);
%     Icolm = Im(:,pickedj);
%     palpha=rs(Icolp)./col(Icolp);
%     maxrcol = max(palpha);
%     malpha = rs(Icolm)./col(Icolm);
%     minrcol = min(malpha);
% 
%     if sum(Icolm)==0
%         inc = maxrcol;
%     elseif sum(Icolp)==0
%         inc = minrcol;
%     elseif maxrcol < minrcol
%         inc = maxrcol;
%     else
% %        active1 = plaha(minrcol < plaha);
% %        active2 = malpha(malpha < maxrcol);
% %        active = sort([active1,active2],'ascend');
% %        AN_r+r(i)*A(i,:)'-active*A(:,i)'*col
%       %  inc = spiecewise(A,b,s*I(:,pickedj),x);
%   %    inc = bisect2(minrcol,maxrcol,rs,col,1e-1);
%              alphSort=sort([minrcol;plaha(minrcol < plaha);malpha(malpha < maxrcol);maxrcol]);
% %        for j = 1:size(alphSort,1)
% %            rj=b-alphSort(j)*col;
% %              rj(rj<0)=0;
% %            df=-col'*rj;
% %            if df >0
% %                inc = 0.5*(alphSort(j)+alphSort(j-1));
% %                 break;
% %            end
% %        end
% %     end
% 
%      Rrepmat = repmat(b,1,size(alphSort,1)) - row*alphSort';
%      Rrepmat(Rrepmat<0)=0;
%      df=-col'*Rrepmat;
%      loc=find(df>0);
%        inc = 0.5*(alphSort(loc)+alphSort(loc-1));
 inc = alpha*( col' * r ) / Acol(pickedj);
%     end


    x(pickedj) = x(pickedj) + inc;
    rst = rs - inc*col;
    r=rst;
   %    if mod(iter,100)==0
       % r( r < 0) = 0;
       r=(r+abs(r))/2;
       I=r>0;
       Ic = I - I_old;
       rs(Ic==-1) = -1*rs(Ic==-1);
       aIc=abs(Ic);
        At_r = At_r + A(aIc,:)'*rs(aIc);
        At_r = At_r-A(I,:)'*col(I);
       rs = rst;
       I_old=I;
 %      At_r = A' * r;
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
function xc=bisect2(a,b,r,row,tol)
fa=f(r,row,a);
fb=f(r,row,b);    
if sign(fa)*sign(fb) >= 0
    error('f(a)f(b)<0 not satisfied!') %ceases execution
end
while (b-a)/2>tol
    c=(a+b)/2;
    fc=f(r,row,c);
    if fc == 0              %c is a solution, done
        break
    end
    if sign(fc)*sign(fa)<0  %a and c make the new interval
        b=c;
        fb=fc;
    else%c and b make the new interval
        a=c;
        fa=fc;
    end
end
xc=(a+b)/2;
end
function fx = f(r,row,x)
    r= r-x*row;
    r(r<0)=0;
    fx = row'*r;
end
end