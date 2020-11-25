% only linear equation
function [xk,fk,dh,rkk]=ssqr4(x0,A,b)
tol=0;
[m,n]=size(A);
r=b-A*x0;
r0=r;
r0(r0<0)=0;
Ar0=norm(A'*r0);
I=find(r>=tol);
AI=A(I,:);
if m==n
    hk=AI\r(I);
else
   % hk=pinv(AI)*r(I);
   [hk,flag] = pcg(AI'*AI,AI'*r(I),1e-10,150);
end
xk=x0+hk;
% ÏÂ½µÁ¿Au
dh=A*hk;
rk=b-A*xk;
rkk=rk;
rkk(rkk<0)=0;
fk=0.5*(rkk'*rkk);
