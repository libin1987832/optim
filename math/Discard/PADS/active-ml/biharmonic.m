function [ x, s, Aopt] = biharmonic( n, repeat);
% biharmonic equation
% input:
%   n        ... initial grid size (e.g. n=10)
%                (number of variables will be n*n)
%   repeat:  ... how often should grid size be doubled (e.g.repeat=3) 
%                (repeat = 0 means solve for given n only)
% output:
% x,s:       ... optimal solution (final grid)
% Aopt:      ... active set at optimum
% call:   [ x, s, Aopt] = biharmonic( n, repeat);
  
tstart = cputime;   % take the time 
Ain = 1:n;          % initial active set
cnt = -1;           % initialize counter

while cnt < repeat;
  disp(['n= ',int2str(n)]);
  [Dn, Fn, Bn] = set_biharm( n);  % generate data for current n, then solve
  [x,s,it,Aopt] = qp_bnd(Dn, -reshape(Fn,n*n,1), reshape(Bn, n*n,1), Ain);
  Ain = sparse(n,n,1);       % initialize new active set
  Ain( Aopt) = ones( length( Aopt),1);
  Ain = kron( Ain, ones(2)); % generate guess on active set at refined grid
  Ain = find( Ain > 0);
  cnt = cnt + 1;
  tcont = cputime;
  disp(['total time (Sec.): ', num2str(tcont-tstart)]);disp(' ');
  n = 2*n;                   % double size
  end;



