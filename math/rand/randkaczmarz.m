function [x0, resvec, xvec, randp]=randkaczmarz(A, b, x0, tol, maxit)
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
prob = Ah' / F;

resvec = zeros(1, maxit);
xvec = zeros(n, maxit);
loopcount = 1;
randp = randsrc(maxit, 1, [number;prob]);
randp = repmat(number, maxit, 1);
r = b - A * x0;
resvec( loopcount )= norm(r);
xvec(:, loopcount )= x0;
while resvec( loopcount ) >= tol && loopcount < maxit
    loopcount = loopcount + 1;
    t = r(randp(loopcount)) / Ah(randp( loopcount ));
    x1 = x0 + t * A(randp(loopcount),:)';
    r = b - A * x1;
    resvec( loopcount )= norm(r);
    xvec(:, loopcount )= x1;
    x0 = x1;
end
resvec = resvec(resvec > eps );
xvec = xvec(:, resvec > eps );

