% A=diag(floor(rand(5,1)*10)) + tril(floor(rand(5,5)),-1)
% inv(A)
m=10;
n=3;
A=rand(m,n);
Ac=diag(sqrt(sum(A.*A,1)));
minAc=min(sum(A.*A,1));
B=A*Ac;
min(eig(A'*A))
minAc*min(eig(B'*B))
