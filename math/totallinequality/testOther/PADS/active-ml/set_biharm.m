function [Dn, Fn, Bn] = set_biharm( n);
% generates discretisation of biharmonic (Dn)
% linear term (Fn) and upper bound (Bn) for
% biharmonic problem as in the paper.  
% call:  [Dn, Fn, Bn] = set_biharm( n);
  
% first coefficient nmatrix Dn  
s = zeros(1,n); s(1) = 20; s(2) = -8; s(3) = 1; 
t = s;
T = toeplitz( s, t);
T = sparse(T);
T(1,1)= 21; T(n,n)=21;
Th = T+speye(n);
S = toeplitz([-8 2 zeros(1,n-2)]);
S = sparse(S);
Dn = kron(speye(n),T);
e = toeplitz([0 1 zeros(1,n-2)]);
e = sparse(e);
Dn= Dn + kron(e,S);
e = toeplitz([0 0 1 zeros(1,n-3)]);
e = sparse( e);
Dn = Dn + kron(e, speye(n));
Dn(1:n,1:n)=Th;
n2=n*n;
Dn(n2-n+1:n2,n2-n+1:n2)=Th;
Dn = (n+1)^4 * Dn;               % scale 

% now force matrix Fn and upper bounds Bn
Fn = zeros(n);
Bn = zeros(n);
for i = 1:n;
for j = 1:n;

  x = i/n; y = j/n;

  Fn( i, j) = -60*(1-x*x)*y*exp( -7*(x-.9)^2 - 4*(y-.1)^2) ...
      + 100*x*(1-y)*exp( -3*(x-.2)^2 - 6*(y-.8)^2);
 
  Bn( i, j) = .00004;
end; end;




