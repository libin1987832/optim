% x0 initial value Q R is decompl A b
% xk out value 
%r0=(b-Ax0)+ rk=(b-Axk)+ 
%fm=b-Axk fr=norm(yk1)^2 
function [xk,r0,rk,fk,fm,fr,l1,l2]=FM(x0,Q,R,A,b)
[m,n]=size(A);
r=b-A*x0;
r0=r;
z0=(-r);
z0(z0<0)=0;
r(r<0)=0;
fr=0.5*(r'*r);
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
xk=x0+uk;

fk=b-A*xk;
fm=fk;
fk(fk<0)=0;
rk=fk;
fk=0.5*fk'*fk;

Lk=-fm-z0;
lfk=0.5*Lk'*Lk;
ukk=A*uk-r;
ufk=0.5*ukk'*ukk;
l1=lfk-fk;
l2=fr-lfk;
