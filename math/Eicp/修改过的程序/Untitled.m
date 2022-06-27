n=5;

c=rand(n,1);
c1=rand(n,1);
c=c-c1;

b=rand(n,1);
b1=rand(n,1);
b=b-b1;

tic
b0=find(c>b)
toc

tic
b1=c>b
toc