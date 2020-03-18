A=rand(4,2)*10;
[q,r]=qr(A);
q=q(:,[1,2]);
Aq=q*q'*diag([1,1,0,0]);
[f,v]=eig(Aq)