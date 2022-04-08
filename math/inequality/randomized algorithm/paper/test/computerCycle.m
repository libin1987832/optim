function skp=computerCycle(A,N,type)
[m,n]=size(A);
NA=N*A;
if type == 1
    ATA=A'*A;
    D1=diag(diag(ATA));
    ATA=A'*N*A;
    D2=diag(diag(ATA));
    U=triu(ATA,1);
    L=tril(ATA,-1);
    DL=D1+lambda*L;
    B=DL*pinv(A);
    skp=eye(m)-lambda*pinv(B)*A'*N;
else
    skp=1;
    for i = 1:n
        Ii=eye(m)-NA(:,i)*NA(:,i)'/sA(i);
        skp=Ii*skp;
    end
end