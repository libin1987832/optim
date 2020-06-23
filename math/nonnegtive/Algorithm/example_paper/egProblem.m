m=1000;
ratio=0.001;
n=ceil(ratio*m);
A=rand(m,n);
xs=rand(n,1);
b=rand(m,1);
r=b-A*xs;
I1=find(r>1e-10);
%b1=pinv(A)*(b(I1)-A(I1,I1)*xs(I1));
b1=A(I1,:)*xs;
b(I1)=b1;
x0=zeros(n,1);
[xk6,fk6,xkArr6,countF6,countN6]=hybrid6(x0,A,b,10);
[xk2,fk2,xkArr2,countF2,countN2]=hybrid2(x0,A,b);
[xk1,fk1,xkArr1,countF1,countN1]=hybrid1(x0,A,b);
[fqf(b,A,xs) fqf(b,A,xk6) fqf(b,A,xk2) fqf(b,A,xk1)]
[norm(xk6-xs) norm(xk2-xs) norm(xk1-xs)]