% x0 initial value Q R is decompl A b
% xk out value 
%r0=(b-Ax0)+ rk=(b-Axk)+ 
%fm=b-Axk fr=norm(yk1)^2 
function [xk,r0,rk,fk,fm,fr]=FM(x0,Q,R,A,b)
[m,n]=size(A);
r=b-A*x0;
r(r<0)=0;
r0=r;
fr=0.5*(r'*r);
% compute min increase
% uk=R\(Q'*r);
Qr=Q'*r;
xk=x0;
for i=n:1
    for j=n:i+1
        xk(i)=Qr(i)/R(i,i);
    end
end
% mm=A*uk-r;
% fm=0.5*mm'*mm;
xk=x0+uk;
fk=b-A*xk;
fm=fk;
fk(fk<0)=0;
rk=fk;
fk=0.5*fk'*fk;
% z0=A*x0-b;
% z0(z0<0)=0;
% u0=A*x0-z0
% 
% zk=A*xk-b;
% zk(zk<0)=0;
% uk=A*xk-zk