% x0 initial value Q R is decompl A b
% xk out value
%r0=(b-Ax0)+ rk=(b-Axk)+
%fm=b-Axk fr=norm(yk1)^2
function [xk,r0,rk,fk,fm,fr]=FMGS(x0,A,b,r)
[m,n]=size(A);
r=(b-A*x0);
r0=r;
r0(r0<0)=0;
fr=0.5*(r'*r);
z0=-r;
z0(z0<0)=0;
bz=b+z0;
xk=x0;
for j=1:3
    for i=1:n
        rk=A*xk-bz;
        c=A(:,i);
        xk(i)=xk(i)-(c'*rk)/(c'*c);
    end
end
fk=b-A*xk;
fm=fk;
fk(fk<0)=0;
rk=fk;
fk=0.5*(fk'*fk);