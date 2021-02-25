clear 
clc
m1 = 5000; n = 1000;
A1 = abs(sprand(m1,n,0.1,1/100));
A2 = rand(m1,n);
%A2(1:floor(m1/3),1:n) = 0;
p =rand(n,1);
t = sparse(p);
A4 = A2 + A2;
% 1
disp 'A1 * A1 * p'
tic 
A1' * (A1 * p);
toc

A3 = A2 + A2;
%2
disp 'A1 * p'
tic 
A1 * p;
toc

A3 = A2+A2;
% 2
disp 'Ap = A1 * p A1 * Ap'
tic
Ap = A1 * p;
A1' * Ap;
toc
A2;
% 3
disp 'A1 * (A1 * t)'
tic 
A1' * (A1 * t);
toc
% 4
A1;
disp 'A2 * p'
tic 
A2 * p;
toc
% 5
disp 'A1 * p'
tic 
A1 * p ;
toc
A2 = A2+A2;
disp 'Ap = A1 * A1  Ap * p'
tic 
Ap = A1' * A1;
Ap * p;
toc
A2 = A2 + A2;
% 3
disp 'A1 * (A1 * p)'
tic 
A1' * (A1 * p);
toc
A2 = A2 + A2;
% 3
disp '(A1 * A1) * p'
tic 
(A1' * A1) * p;
toc