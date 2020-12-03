function [x0, resvec, randp]=randkaczmarz(A, b, x0, tol, maxit)
[m,n]=size(A);

if nargin < 4   
    maxit = 1e4;
end


if size(b, 1) ~= m || size(b, 2) ~= 1
    error('The size of A and b do not match')
elseif size(x0, 1) ~= n || size(x0, 2) ~= 1
    error('The size of x0 does not match the problem')
end

number = 1 : m;
As = A .* A;
Ah = sum(As, 2);
F = sum(Ah);
prob = Ah / F;

resvec = zeros(1, maxit);
loopcount = 0;
randp = randsrc(maxit, 1, [number;prob]);

while resvec( loopcount ) >= tol && loopcount < maxit
    loopcount = loopcount + 1;
    r = A * x0 - b;
    resvec( loopcount )= norm(r);
    t = - r(randp(loopcount)) / Ah(randp( loopcount ));
    x0 = x0 + t * A(randp(loopcount),:)';
end
resvec(resvec < eps ) = 0;

