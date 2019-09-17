function [xk,fk,y]=ssqr2(x0,A,b)
tol=0;
[m,n]=size(A);
r=b-A*x0;
I=find(r>=tol);
AI=A(I,:);
AII=AI'*AI;
bII=AI'*(r(I));
[U,S,V]=svd(AII);
rnkd=length(find(diag(S) >= 1e-14));
Udd=U(:,1:rnkd);
Sdd=S(1:rnkd,1:rnkd);
Vdd=V(:,1:rnkd);
hk=Vdd*(Sdd\(Udd'*bII));

xk=x0+hk;
dh=A*hk;
rk=b-A*xk;
I1=find(rk>=tol);
y=isequal(I,I1);

if ~y && rnkd<n
    ai=r./dh;
    aa=min(ai(ai>tol));
    xk=x0+aa*hk;
    rk=b-A*xk;
    I1=find(rk>=tol);
    y=all(ismember(I,I1));
end

rk=rk(I1);
fk=0.5*rk'*rk;
   






























































































































































































