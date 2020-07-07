% x0 initial value Q R is decompl A b
% xk out value 
%r0=(b-Ax0) rk=(b-Axk)  
function [uk,r0,rk]=FixedM(x0,Q,R,A,b)
[m,n]=size(A);
r=b-A*x0;
r0=r;
r(r<0)=0;
Qr=Q'*r;
% compute min increase
uk=R\Qr;
% uk=x0;
% for i=n:-1:1
%     for j=n:-1:(i+1)
%         Qr(i)=Qr(i)-uk(j)*R(i,j);
%     end
%     uk(i)=Qr(i)/R(i,i);
% end
% mm=A*uk-r;
% fm=0.5*mm'*mm;
rk=r-A*uk;

