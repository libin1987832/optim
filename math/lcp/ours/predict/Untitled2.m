n=3;
% A=rand(n);
% xs=rand(n,1);
% save('tt1w','A','xs');

 load('tt','A','xs')
q=A*xs;
A\q
xs
D = diag(diag(A));%求A的对角矩阵
L = tril(A,-1);%求A的下三角矩阵,不带对角线
U = triu(A,1);%求A的上三角矩阵
DLU=-inv(D+L)*U;


[x,n,xa]=guaseidel(A,q,zeros(n,1),1e-5,10);
xa;
x;
xad=xa(:,1:8)-xa(:,2:9);
xro=xad(:,2:8)./xad(:,1:7);
abs(xro);
l=max(eig(DLU));
xp=xa(:,4)+(xa(:,5)-xa(:,4))/(1-l)
xpp=predict2(xa(:,4),xa(:,5),xa(:,6))
