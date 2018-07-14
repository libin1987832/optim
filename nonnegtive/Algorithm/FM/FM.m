function [xk,fk]=FM(x0,Q,R,A,b)
r=b-A*x0;
xk=R\(Q\r);
fk=b-A*x0;
fk(fk<0)=0;