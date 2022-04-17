m=5;n=3;
A=rand(m,n)
[u,v,s]=svd(A);
u(:,1:n)*v(1:n,:)*s'
Ai=s*diag(1./diag(v))*u(:,1:n)'
N=diag([0,1,0,1,0]);
Ai*N*A
[u1,v1,s1]=svd(Ai*N*A);
ATA=A'*A;
ATAi=inv(ATA)
ATABA=ATAi*A'*N*A
ATABA2=Ai*N*A
[V,D]=eig(ATABA)
eig(V)
% V*D
% ATABA*V
% V'*V
% [u1,v1,s1]=svd(ATABA)
% v1