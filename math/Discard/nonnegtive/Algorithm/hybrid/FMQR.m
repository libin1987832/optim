% x0 initial value Q R is decompl A b
% xk out value 
%r0=(b-Ax0)+ rk=(b-Axk)+ 
%fm=0.5*norm(A*uk-r)^2 fr=norm(yk1)^2 
function [xk,fk]=FM1(x0,Q,R,A,b)
r=b-A*x0;
r(r<0)=0;
% compute min increase
uk=R\(Q'*r);

xk=x0+uk;
rk=b-A*xk;
rk(rk<0)=0;

fk=0.5*rk'*rk;

% z0=A*x0-b;
% z0(z0<0)=0;
% u0=A*x0-z0
% 
% zk=A*xk-b;
% zk(zk<0)=0;
% uk=A*xk-zk