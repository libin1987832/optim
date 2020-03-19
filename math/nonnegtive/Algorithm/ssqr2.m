% Pinar method
function [xk,fk,y]=ssqr2(x0,A,b)
% Pinar method
tol=0;
[m,n]=size(A);
r=b-A*x0;
r0=r;
r0(r0<0)=0;
I=find(r>=tol);
%提取子矩阵判断是否正定
AI=A(I,:);
% AII=AI'*AI;
% bII=AI'*(r(I));
% sparse for svds densy for svd
% [U,S,V]=svds(AII);
% %找大于0的所有值 判断正定
% rnkd=length(find(diag(S) >= 1e-20));
%
% Udd=U(:,1:rnkd);
% Sdd=S(1:rnkd,1:rnkd);
% Vdd=V(:,1:rnkd);
% hk=Vdd*(Sdd\(Udd'*bII));
hk=AI\r(I);
xk=x0+hk;
% 下降量Au
dh=A*hk;
rk=b-A*xk;
rkk=rk;
rkk(rkk<0)=0;
% 牛顿法的结果是否保持积极集不变
I1=find(rk>=tol);
y=isequal(I,I1);
% 寻找最小的那个改变的值 -a'(y+ai*h)+c=0 (-Ay+c)(r)+ai*(dh)
% if ~y && rnkd<n
if ~y
    d = eig(AI'*AI);
    tol2=length(d)*eps(max(d));
    isposdef = all(d) > tol2;
    issemidef = all(d) > -tol2;
end

%if ~y && rnkd<n
if ~y && isposdef
    % nonsingular with exact arithtic (there is some issue)
    %         ai=r./dh;
    %         aa=min(ai(ai>tol));
    %         rt=r-aa*dh;
    %         a=r(rt>0)'*dh/(dh'*dh);
    aa=bisect(r,dh);
    % aa=plusoperator(A,b,dh,x0);
    xk=x0+aa*hk;
    rk=b-A*xk;
    rkk=rk;
    rkk(rkk<0)=0;
    I1=find(rk>=tol);
    % I 中的元素是否在 I1中 如果在其中就是true 如果没有就是false 如果全是true则为1 否则全为0
    y=all(ismember(I,I1));
    
    %else if ~y && rnkd==n
else if ~y && issemidef
        %  singular first breakpoints
        ai=r./dh;
        aa=min(ai(ai>tol));
        xk=x0+aa*hk;
        rk=b-A*xk;
        rkk=rk;
        rkk(rkk<0)=0;
        I1=find(rk>=tol);
        % I 中的元素是否在 I1中 如果在其中就是true 如果没有就是false 如果全是true则为1 否则全为0
        y=all(ismember(I,I1));
    end
end
% 截取违背的不等式
rk=rk(I1);
% 返回目标函数的值
fk=0.5*rk'*rk;
