% x0 initial value Q R is decompl A b sparse
% xk out value 
%r0=(b-Ax0)+ rk=(b-Axk)+ 
%fm=b-Axk fr=norm(yk1)^2 
function [xk,r0,rk,fk,fm,fr]=FMS(x0,Q,R,A,b)
[m,n]=size(A);
r=b-A*x0;
r(r<0)=0;
r0=r;
fr=0.5*(r'*r);
Qr=Q'*r;

% compute min increase
%uk=R\Qr;
uk=x0;
for i=n:-1:1
%     [ix,iy,iv]=find(R);
    for j=n:-1:(i+1)
        Qr(i)=Qr(i)-uk(j)*R(i,j);
    end
    uk(i)=Qr(i)/R(i,i);
end
% mm=A*uk-r;
% fm=0.5*mm'*mm;
xk=x0+uk;
fk=b-A*xk;
fm=fk;
fk(fk<0)=0;
rk=fk;
fk=0.5*fk'*fk;
