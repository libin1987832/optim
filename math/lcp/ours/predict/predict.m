function [x]=predict(A,x0,q,iter)
[m,n]=size(A);
b=-1*computAC(A,q);
rq=b;
for i=1:iter
rq=computDLU(A,rq);
rq=rq+b;
end
x=rq;