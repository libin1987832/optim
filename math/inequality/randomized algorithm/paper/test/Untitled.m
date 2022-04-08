A=rand(5,3);
b=rand(5,1);
x00=rand(3,1);
x1=x00;
N=diag([1,0,1,0,1]);
lambda = 1.1;
rx0=b-A*x00;
for i=1:3
    x1(i)=x1(i)+lambda*A(:,i)'*N*(b-A*x1)/(A(:,i)'*A(:,i));
    rx1=b-A*x1
    rx1t=(eye(5)-lambda*(N*A(:,i))*(N*A(:,i))'/(A(:,i)'*A(:,i)))*rx0
end
ATA=A'*A;
D1=diag(diag(ATA));
ATA=A'*N*A;
D2=diag(diag(ATA));
U=triu(ATA,1);
L=tril(ATA,-1);
DL=D1+lambda*L;
y00=x00;
x01=-inv(DL)*(lambda*U+lambda*D2-D1)*x00+lambda*inv(DL)*A'*N*b;
x02=-inv(DL)*(lambda*U+lambda*D2-D1)*x01+lambda*inv(DL)*A'*N*b;
r01=b-A*x01;
r00=b-A*x00;
r02=b-A*x02
B=DL*pinv(A);
C=eye(5)-lambda*pinv(B)*A'*N;
r11=C*r00;
r12=C*r11
%  lambda*A'*N*r00
%  [U,S,V]=svd((N*A)');
%  n1=(N*r00)'*V(1,:)'
%   n2=(N*r01)'*V(1,:)'
%    n3= (N*r02)'*V(1,:)'
%    n2/n1
%    n3/n2
%  (D1+lambda*L)*(x00-x1)
% [U,S,V]=svd(A);
% U=U(:,1:3);
% S=S(1:3,1:3);
% lambda*V*S*U'*N*r00
% (D1+lambda*L)*V*diag([1/S(1,1),1/S(2,2),1/S(3,3)])*U'*(r01-r00)

% x01
% x1
% r01=b-A*x01;
% r00=b-A*x00;
% 2*lambda*A'*N*r01
% lambda*A'*N*r00-(D1+lambda*L)*(x00-x1)



% for i=1:200
% y01=-inv(DL)*(lambda*U+lambda*D2-D1)*y00+lambda*inv(DL)*A'*N*b;
% y00=y01;
% end
% x01=-inv(DL)*(lambda*U+lambda*D2-D1)*x00+lambda*inv(DL)*A'*N*b;
% x01
% y00
% x02=-inv(DL)*(lambda*U+lambda*D2-D1)*x01+lambda*inv(DL)*A'*N*b;
% [U,S,V]=svd((N*A)');
% r1=b-A*x00;
% r2=b-A*x01;
% r3=b-A*x02;
% r1p=N*(b-A*x00);
% r2p=N*(b-A*x01);
% r3p=N*(b-A*x02);
% i1=r1'*V(:,1);
% i2=r2'*V(:,1);
% i3=r3'*V(:,1);
% i1p=r1p'*V(:,1);
% i2p=r2p'*V(:,1);
% i3p=r3p'*V(:,1);
% isp=(b-A*y01);
% isp'*V(:,1)
% norm(-inv(DL)*(lambda*U+lambda*D2-D1))
% [i1,i2,i3,i1/i2,i2/i3;i1p,i2p,i3p,i1p/i2p,i2p/i3p]


