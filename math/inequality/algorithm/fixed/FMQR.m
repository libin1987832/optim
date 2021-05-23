% x0 initial value Q R is decompl A b
% xk out value 
%r0=(b-Ax0)+ rk=(b-Axk)+ 
%fm=0.5*norm(A*uk-r)^2 fr=norm(yk1)^2 
function [xk,rk]=FMQR(x0,Q,R,A,b,rk)
rk(rk<0)=0;
% compute min increase
QTrk=(Q'*rk);

uk=R\QTrk;
xk=x0+uk;
rk=b-A*xk;
