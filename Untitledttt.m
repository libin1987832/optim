a=1;
for i=1:10
a1=0.5*(1+sqrt(4*a^2+1));
(a-1)/a1
a=a1;
end

% m=5;n=3;
% A=rand(m,n)
% [u,v,s]=svd(A);
% u(:,1:n)*v(1:n,:)*s'
% Ai=s*diag(1./diag(v))*u(:,1:n)'
% N=diag([0,1,0,1,0]);
% Ai*N*A
% [u1,v1,s1]=svd(Ai*N*A);
% ATA=A'*A;
% ATAi=inv(ATA)
% ATABA=ATAi*A'*N*A
% ATABA2=Ai*N*A
% [V,D]=eig(ATABA)
% <<<<<<< HEAD
% V*D*pinv(V)
% =======
% eig(ATABA2)
% svd(ATABA2)
% eig(u(:,1:n)'*N*u(:,1:n))
% eig(u(:,1:n)'*(eye(m)-N)*u(:,1:n))
% >>>>>>> bcfab28382feef822e224ef4b6a0679c77ddad0f
% V*D
% ATABA*V
% V'*V
% [u1,v1,s1]=svd(ATABA)
% v1