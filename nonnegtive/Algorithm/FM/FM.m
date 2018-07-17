function [xk,fk]=FM(x0,Q,R,A,b)
r=b-A*x0;
r(r<0)=0;
uk=R\(Q\r);
xk=x0+uk;
fk=b-A*xk;
fk(fk<0)=0;
fk=0.5*fk'*fk;