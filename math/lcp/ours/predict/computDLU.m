function xk=computDLU(A,x0)
[m,n]=size(A);
xk=zeros(n,1);
for j=1:n
    xk(j)=x0(j)-(A(j,:)*[xk(1:(j-1));x0(j:n)])/A(j,j );
end
%