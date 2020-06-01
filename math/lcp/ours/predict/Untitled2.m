n=3;
A=rand(n);
xs=rand(n,1);
q=A*xs;
A\q
xs
D = diag(diag(A));%求A的对角矩阵
L = tril(A,-1);%求A的下三角矩阵,不带对角线
U = triu(A,1);%求A的上三角矩阵
DLU=-inv(D+L)*U;
eig(DLU)
[x,n,xa]=guaseidel(A,q,zeros(n,1),1e-5,10);
xa;
x;