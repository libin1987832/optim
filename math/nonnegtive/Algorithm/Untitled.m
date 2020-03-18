% A=rand(4,2)*10;
% [q,r]=qr(A);
% q=q(:,[1,2]);
% Aq=q*q'*diag([1,1,0,0])
% [f,v]=eig(Aq)
A=rand(4,2)*10;
As=sym(floor(A));
[q,r]=qr(As);
q=q(:,[1,2]);
Aq=q*q'*diag([1,1,0,0]);
Aq;
[v,j]=jordan(B);
lj=logical(vpa(j>0.99991));
jj=diag([lj(1,1),lj(2,2),lj(3,3),lj(4,4)]);
vpa(v*jj*inv(v))
% [f,v]=eig(Aq)