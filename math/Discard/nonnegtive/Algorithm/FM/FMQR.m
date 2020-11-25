% x0 initial value Q R is decompl A b
% xk out value 
%r0=(b-Ax0)+ rk=(b-Axk)+ 
%fm=0.5*norm(A*uk-r)^2 fr=norm(yk1)^2 
function [xk,r0,rk,z,fk,fr,fu,fz,fd]=FMTD(x0,Q,R,A,b)
r=b-A*x0;
r(r<0)=0;
r0=r;
fr=0.5*r'*r;
% compute min increase
uk=R\(Q'*r);
Au=A*uk;
fu=0.5*Au'*Au;

xk=x0+uk;
fk=b-A*xk;
fk(fk<0)=0;
rk=fk;
fk=0.5*fk'*fk;

z=r-A*uk;
fz=0.5*z'*z;
fd=fz-fk;
% z0=A*x0-b;
% z0(z0<0)=0;
% u0=A*x0-z0
% 
% zk=A*xk-b;
% zk(zk<0)=0;
% uk=A*xk-zk