function [ x, s, iter, Aopt] = qp_bnd( Q, d, b, A);

% solves: min d'x + (1/2) x'Qx  subject to: x  <= b
% Q is assumed to be symmetric positive definite
% KKT:    Qx + s + d = 0;  x <= b, s >=0, s'(x-b) = 0
% input:  (Q,d,b) problem data,
%         A ... guess on initial active set, e.g. A=(1:n)�; 
% output: (x,s) optimal solution, iter= # iterations
%         Aopt = active set at opt. sol. 
% call:   [ x, s, iter, Aopt] = qp_bnd( Q, d, b, A);

Q = sparse(Q);  
n = length( d);                  % problem size
k = 0;                           % iteration count
done = 0;                        % not yet done
tol =1e-12;                      % tolerance for feasibility

disp(' iter inact '  );          % display some info on screen 

% main loop 
while done < 1;                  % while not optimal
     k = k + 1;                  % start a new iteration
% solve system KKT(A):
     I = ones(n,1);              % compute I, the complement of A
     I(A) = zeros( length(A),1); 
     I = find(I>0);              % complement of A
     x = b;                      % x( A) = b( A); 
     s = zeros( n,1);            % s( I) = 0
     if length( I) > 0;          
%         ltri = sparse(chol(full(Q( I,I))));  % cholesky factor
         id = I;
         p = symrcm(Q(I,I));
	 [dummy,p_inv]=sort(p);
	 I1 = id(p);
         ltri = sparse(chol((Q( I1,I1))));  % cholesky factor
         rhs = -d( I1);
         if length(A) > 0;       % update right hand side rhs
           rhs = rhs - Q( I1,A)*x( A);   
         end; 	   
         xtmp = (ltri') \ rhs; 
         x( I) = ltri \  xtmp   % solve for inactive variables	   
	 x( I) = x( id(p_inv));
     end;                               
     if length( A)>0;            % backsubsitute for s(A), if |A|>0 
         s( A) = -d( A) - Q( A, :)* x  
     end;

     A = [find(x>b); find(s>0)]; % active constraints of new solution

     done = (max(x-b)<= tol) & (min(s)>= -tol); 

     fprintf(' %3.0d %6.0f \n',[k  (n-length(A)) ]);
     if k > 100; done = 1;       % emergency exit to avoid cycling
         disp(' max number of iterations reached.');
     end;
end;                 % end while
iter = k; 
Aopt = A;












