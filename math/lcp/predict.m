function [x]=predict(A,x0,q)
[m,n]=size(A);
b=computAC(A,q);
for i=1:m
    for j=1:i-1
        bd=bd+A(i,j)*b(j);
    end
    b(i)=(q(i)-bd)/A(i,i);
end
C=