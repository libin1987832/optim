A1=abs(sprand(m1,n,0.1,1/100));
A2 = rand(m1,n);
p =rand(n,1);
t = sparse(p);
tic
A1 * p;
toc
tic 
A2 * p;
toc