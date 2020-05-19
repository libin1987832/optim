function b=computAC(A,c)
[m,n]=size(A);
b=zeros(n,1);
for i=1:m
    bd=A(i,1:i-1)*b(1:i-1);
    b(i)=(c(i)-bd)/A(i,i);
end
