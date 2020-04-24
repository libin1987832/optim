function [xk,rk,fk,f0,lambe]=ssqr(x0,A,b)
% take the binary search step length
r=b-A*x0;
I=(r>=0);
AI=A(I,:);

%[Q,R]=qr(AI);
%h=r(I);
%u=R\Q'*h;

AII=AI'*AI;
bII=AI'*(r(I));
[U,S,V]=svds(AII);

rnkd=length(find(diag(S) >= 1e-10));
Udd=U(:,1:rnkd);
Sdd=S(1:rnkd,1:rnkd);
Vdd=V(:,1:rnkd);
u=Vdd*(Sdd\(Udd'*bII));
% 内置算法 有时候会比SVD要快
% u=AI\r(I);

r(r<0)=0;
f0=0.5*(r'*r);

lambe=1;
xk=x0+lambe*u;
r=b-A*xk;
r(r<0)=0;
fk=0.5*(r'*r);

while fk>f0
    lambe=0.5*lambe;
    xk=x0+lambe*u;
    r=b-A*xk;
    r(r<0)=0;
    fk=0.5*r'*r;
end
rk=r;




