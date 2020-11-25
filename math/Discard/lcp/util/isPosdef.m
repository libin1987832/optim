function [symmetric,opsdef]=isPosdef(M)
s=issparse(M);
symmetric = issymmetric(M);
if s
d=eigs(M);
else
d = eig(M);
end
tol=length(d)*eps(max(d));
opsdef = all(d > tol);
if opsdef==0 && all(d > -tol)
  opsdef=2;
end