function skp=computerMSE(A,N,c)
[m,n]=size(A);
sA=sum(A.*A);
NA=N*A;
snA=sum(NA.*NA);
ssA=sum(sA);
skp=0;
for i = 1:n
        pa=eye(m)-NA(:,i)*NA(:,i)'/sA(i);
    if c==1
       kp=kron(pa,pa)/n;
    else
       kp=sA(i)*kron(pa,pa)/ssA;
    end
     skp=skp+kp;
end