A1 = abs(sprand(m1,n,0.1,1/100));
A2 = rand(m1,n);
%A2(1:floor(m1/3),1:n) = 0;
p =rand(n,1);
t = sparse(p);
A4 = A2 + A2;
tic 
A1' * (A1 * p);
toc
A3 = A2+1;
tic
Ap = A1 * p;
A1' * Ap;
toc
A2;
tic 
A1' * (A1 * t);
toc
A1;
tic 
A2 * p;
toc
tic 
A1 * p ;
toc