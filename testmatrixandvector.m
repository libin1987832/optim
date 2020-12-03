A = rand(12000,400);
B = rand(400,12000);
% A.*A, 2
f = @() A.^2;
timeit(f)

g =@() norm(A,'fro')^2;
timeit(g)
norm(A,'fro')^2