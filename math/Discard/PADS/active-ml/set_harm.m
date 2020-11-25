function [Q, d, b] = set_harm( k);
% generates symmetric block tridiagonal matrix of order k*k
% with (i,i)-block = toeplitz( 4 -1 0 0 ...), and (i,i+1)-block=-I
% and input for qp  
% call [Q,d,b] = set_harm( k):
  
  a1 = toeplitz([ 4 -1 zeros(1,k-2)]);
  a1 = sparse(a1);
  a2 = toeplitz([ 0 1 zeros(1,k-2)]);
  a2 = sparse( a2);
  Q= kron( speye(k), a1)+ kron( a2, -speye(k));

Bn = zeros(k);
for i = 1:k;
for j = 1:k;
  
  x = i/(k+1); y = j/(k+1);
  Bn( i, j) = sin( 3.2*x) * sin(3.3*y);
  end; end;

b = reshape( Bn, k*k, 1);
d = -ones(k*k,1);