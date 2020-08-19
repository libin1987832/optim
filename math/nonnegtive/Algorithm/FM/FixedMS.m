% x0 initial value Q R is decompl A b
% xk out value 
%r0=(b-Ax0) rk=(b-Axk)  
function [xk,rk]=FixedMS(x0,Q,U,A,b,r)
[m,n]=size(A);
r(r<0)=0;
Qr=Q'*r;
% compute min increase
%uk=R\Qr;
 
% uk=x0;
% uk=zeros(n,1);
% for i=n:-1:1
%      uk(i)=(Qr(i)-R(i,i+1:end)*uk(i+1:end))/R(i,i);
% end



uk(1,n)=Qr(end)./U(n,n);
%zeros(1,m);
%bacward substitution
for k=n-1:-1:1 
        x1=1/U(k,k).*(Qr(k)-sum(U(k,k+1:end).*uk(k+1:end)));
        uk(k)=x1;
end
uk=uk';

xk=x0+uk;
rk=b-A*xk;
% mm=A*uk-r;
% fm=0.5*mm'*mm;
%rk=r0-A*uk;

