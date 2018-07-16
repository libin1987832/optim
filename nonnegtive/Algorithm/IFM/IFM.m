function [xk,fk]=IFM(x0,A,b)
r1=b-A*x0;
y1=r1;
y1(y1<0)=0;
[uk,fk]=kyrlov(A,y1);
xk=x0+uk;