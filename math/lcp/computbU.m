function b=computbU(U,b)
[m,n]=size(U);
ub=zeros(m,1)
for i=1:m
    ub(i)=U(i,i:n)*b(i:n);
end