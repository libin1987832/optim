function [xk,fk,y,isP]=ssqr3(x0,A,b)
tol=0;
[m,n]=size(A);
r=b-A*x0;
I=find(r>=tol);
%提取子矩阵判断是否正定
AI=A(I,:);
AII=AI'*AI;
bII=AI'*(r(I));
[U,S,V]=svd(AII);
%找大于0的所有值 判断正定
rnkd=length(find(diag(S) >= 1e-14));
if rnkd == length(I)
    isP=1;
    Udd=U(:,1:rnkd);
    Sdd=S(1:rnkd,1:rnkd);
    Vdd=V(:,1:rnkd);
    hk=Vdd*(Sdd\(Udd'*bII));
    
    xk=x0+hk;
    rk=b-A*xk;
else
    isP=0;
    xk=x0;
end
rk=b-A*xk;
I1=find(rk>=tol);
y=isequal(I,I1);%face is diff,0 diff 1 same
rk=rk(I1);
fk=0.5*(rk'*rk);



