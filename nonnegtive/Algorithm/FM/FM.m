function [xk,fk]=FM(x0,Q,R,A,b)
r=b-A*x0;
r(r<0)=0;
uk=R\(Q\(A'*r));
xk=x0+uk;
fk=b-A*xk;
fk(fk<0)=0;
fk=0.5*fk'*fk;

z0=A*x0-b;
z0(z0<0)=0;
u0=A*x0-z0

zk=A*xk-b;
zk(zk<0)=0;
uk=A*xk-zk