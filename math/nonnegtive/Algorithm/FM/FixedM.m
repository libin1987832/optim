% x0 initial value Q R is decompl A b
% xk out value 
%r0=(b-Ax0) rk=(b-Axk)  
function [xk,rk]=FixedM(x0,Q,R,A,b,r)
[m,n]=size(A);
r(r<0)=0;
Qr=Q'*r;
% compute min increase
uk=R\Qr;
xk=x0+uk;
rk=b-A*xk;
% uk=x0;
% for i=n:-1:1
%     for j=n:-1:(i+1)
%         Qr(i)=Qr(i)-uk(j)*R(i,j);
%     end
%     uk(i)=Qr(i)/R(i,i);
% end
% mm=A*uk-r;
% fm=0.5*mm'*mm;
%rk=r0-A*uk;

