function [xk,zk,fk]=PC(x1,z1,A,b)
e1=A'*(A*x1-z1-b);
r=A*x1-b;
r(r<0)=0;
e2=z1-r;
sum=e1'*e1+e2'*e2;
ae=A*e1-e2;
p=(sum)/(sum+ae'*ae);
xk=x1-p*e1;
zk=z1-p*e2;
fk=b-A*x0;
fk(fk<0)=0;