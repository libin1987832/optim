n=3;
A=rand(n)*1000;
xs=rand(n,1)*10;
q=A*xs;
A\q
xs
[x,n,xa]=guaseidel(A,q,ones(n,1)*10,1.0e-5,50);
x