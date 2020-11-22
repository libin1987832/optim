A = [1 2;3 4];
b = [4;5];
x0 = [0;0];
   
D = A.*A;
D = sum(D,1);

rk0 = b - A * x0;
rk = rk0;
[xk1,rk] = FixedGS(x0,A,b,D,rk,1);

ATA = A'*A;
D = diag(diag(ATA));
L = tril(ATA) - D;
U = triu(ATA) - D;
iDL = inv(D+L);
zk = -rk0;
zk(zk<0) = 0;
ATb = A' * (b+zk);
xk2 = -iDL * ( U * x0 -  ATb);
[xk1 xk2]
