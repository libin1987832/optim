A=rand(3,3);
b=rand(3,1);
x00=rand(3,1);
x1=x00;
N=diag([1,0,1]);
lambda = 1.1;
for i=1:3
    x1(i)=x1(i)+lambda*A(:,i)'*N*(b-A*x1)/(A(:,i)'*A(:,i));
end
ATA=A'*A
D1=diag(diag(ATA));
ATA=A'*N*A
D2=diag(diag(ATA));
U=triu(ATA,1)
L=tril(ATA,-1)
DL=D1+lambda*L;
x01=-inv(DL)*(lambda*U+lambda*D2-D1)*x00+lambda*inv(DL)*A'*N*b;

x1