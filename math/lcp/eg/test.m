n=10000;
A=sprandsym(n,0.4);
x=sprandn(n,1,0.3);
tic 
a1=A(4,:)*x;
f1=toc;
tic
a2=full(A(4,:))*full(x);
f2=toc;
tic
I1=find(A(4,:)>0) & find(x>0);
a3=full(A(4,I1))*full(x(I1));
f3=toc;
fprintf('ff%8.4f,%8.4f,%8.4f\n',f1,f2,f3);