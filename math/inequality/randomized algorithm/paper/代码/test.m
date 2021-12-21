r0=b-A*x0;
rp0=r0;
rp0(rp0<0)=0;
norm_a = A(:,1)'*A(:,1);
inc = 2*A(:,1)'*rp0/norm_a;
r1=r0-inc*A(:,1);
rp1=r1;
rp1(rp1<0)=0;
er1 = error_GS1(:,2);
norm(rp1-er1)
[m,n]=size(A);
N0=diag(r0>0);
det1=(N0*A(:,1)*(N0*A(:,1))')/norm_a;
det2=inc;
Nr1=r0-2*det1*rp0-inc*(eye(m)-N0)*A(:,1);
Nr1(Nr1<0)=0;
norm(Nr1-er1)
AA=A(:,1)*A(:,1)';
AAr1=r0-2*AA/norm_a*rp0;
AAr1(AAr1<0)=0;
norm(AAr1-er1)
% norm(inc*(eye(m)-N0)*A(:,1))