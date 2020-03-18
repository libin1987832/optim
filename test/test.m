% A=qr(rand(3,3));
% [q,r]=qr(A)
% T=q'*diag([1/2,1/3,1])*q;
% T'*T
A=sym([0.5,0.3;1,0.7;1,3;0.3,0.6])
[q,r]=qr(A)
Aq=q(:,[1,2]);
B=diag([1,1,1,1])-Aq*Aq'*diag([1,1,0,0])
[v,j]=jordan(B);
vpa(j,4)
tt=v*diag([0,0,1,1])*inv(v)
vpa(tt,4)
tt=v*inv(v)
vpa(tt,4)