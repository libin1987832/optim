function [x]=predict(A,x0,q)
[m,n]=size(A);
b=computAC(A,q);
rq=b;
for i=1:3
rq=computDLU(A,rq);
rq=rq-b;
end
x=rq;