% x0 initial value Q R is decompl A b
% xk out value 
%r0=(b-Ax0) rk=(b-Axk)  
function [xk,rk]=sfixedM(x0,Q,R,A,b,r)
[m,n]=size(A);
r(r<0)=0;
Qr=Q'*r;
% compute min increase
uk=R\Qr;
xk=x0+uk;
rk=b-A*xk;


