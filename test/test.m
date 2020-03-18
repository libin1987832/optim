% A=qr(rand(3,3));
% [q,r]=qr(A)
% T=q'*diag([1/2,1/3,1])*q;
% T'*T
A=sym([0.5,0.3;1,0.7])
[q,r]=qr(A)