A=rand(3,3);
b=rand(3,1);
x00=rand(3,1);
x1=x00;
N=diag([1,0,1]);
for i=1:3
    x1(i)=x1(i)+A(:,i)'*N*(b-A*x1)/(A(:,i)'*A(:,i));
end
ATA=A'*N*A;
D=diag(diag(ATA));
U=triu(ATA,1);
L=tril(ATA,-1);
DL=D+L;
x01=-inv(DL)*U*x00+inv(DL)*A'*N*b
x1