% Pinar method
function [xk,fk,y]=ssqr2(x0,A,b)
% Pinar method
tol=0;
[m,n]=size(A);
r=b-A*x0;
I=find(r>=tol);
%��ȡ�Ӿ����ж��Ƿ�����
AI=A(I,:);
AII=AI'*AI;
bII=AI'*(r(I));
[U,S,V]=svd(AII);
%�Ҵ���0������ֵ �ж�����
rnkd=length(find(diag(S) >= 1e-14));

Udd=U(:,1:rnkd);
Sdd=S(1:rnkd,1:rnkd);
Vdd=V(:,1:rnkd);
hk=Vdd*(Sdd\(Udd'*bII));

xk=x0+hk;
% �½���Au
dh=A*hk;
rk=b-A*xk;
% ţ�ٷ��Ľ���Ƿ񱣳ֻ���������
I1=find(rk>=tol);
y=isequal(I,I1);
% Ѱ����С���Ǹ��ı��ֵ -a'(y+ai*h)+c=0 (-Ay+c)(r)+ai*(dh)
% if ~y && rnkd<n
if ~y    
    ai=r./dh;
    aa=min(ai(ai>tol));
    xk=x0+aa*hk;
    rk=b-A*xk;
    I1=find(rk>=tol);
    % I �е�Ԫ���Ƿ��� I1�� ��������о���true ���û�о���false ���ȫ��true��Ϊ1 ����ȫΪ0
    y=all(ismember(I,I1));
end
% ��ȡΥ���Ĳ���ʽ
rk=rk(I1);
% ����Ŀ�꺯����ֵ
fk=0.5*rk'*rk;
   