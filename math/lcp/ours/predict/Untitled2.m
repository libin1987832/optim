n=3;
A=rand(n);
xs=rand(n,1);
q=A*xs;
A\q
xs
D = diag(diag(A));%��A�ĶԽǾ���
L = tril(A,-1);%��A�������Ǿ���,�����Խ���
U = triu(A,1);%��A�������Ǿ���
DLU=-inv(D+L)*U;
eig(DLU)
[x,n,xa]=guaseidel(A,q,zeros(n,1),1e-5,10);
xa;
x;