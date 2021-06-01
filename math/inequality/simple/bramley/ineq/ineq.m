%
%   [x,fval,gradient,err,J,D,iters] = ineq(A,x,b,tolerances,dirmeth,linemeth,maxits,printlevel,neq);
%
% Implement Han's method for solving linear inequalities
%
%   Randall Bramley and Beata Winnicka
%   Dept of Computer Science
%   Indiana University - Bloomington
%   24 August 1993
%   Revised 1 October 1993
%   Revised 15 November 1993
%
%   On entry, 
%      x = initial estimate of solution
%      b = rhs
%      A = matrix; all rows are assumed to be for inequalities in this code.
%      tolerances = a two entry array, currently containing the tolerances on
%                   the function/gradient norms and the line search, resp.
%      dirmeth = integer giving method for search direction finding.
%      linemeth = integer giving method for line search.
%      maxits = maximum allowed iterations.
%      printlevel = amount of printing to perform.  Not fully implemented yet.
%                 = 0 perform limited printing (mainly error checks)
%                 = 1 print out function and gradient data from each step
%                 = 2 print out line search and search direction data, and iteration info
%      neq = number of rows of A corresponding to equality conditions.  Note that
%            if all rows of A are for equalities, then one iteration should work but 
%            because of rounding errors more may be taken.
%            *** the equality rows need to be listed LAST in the rows of A***
%
%   On exit,
%      x = solution
%      fval = norm of function value
%      gradient = norm of gradient
%      err = error code; currently only is nonzero if d is not a descent direction.
%      J = array containing active indices from each iteration; see plotindex.m
%      D = array containing search directions from each iteration. (Not used currently)
%      iters = total number of iterations required
%
    function [x,fval,gradient,err,J,iters] = ineq(A,x,b,tolerances,dirmeth,linemeth,maxits,printlevel,neq);
%
%------------------------------------------------------------------------------
% Begin: find initial active index set, residual, function value and gradient.
%------------------------------------------------------------------------------
%
    [m,n] = size(A);
% Number of inequality rows of A:
    ninq = m - neq;
% Index sets corresponding to inequality and equality rows of A:
    INQ = 1:ninq;
    EQ = ninq+1:m;
    INQ = INQ';
    EQ = EQ';
% Note that actually iter should start out at 0, but it is used for indexing
% into J and D, and so must be positive.
    iter = 1;
    r = b - A*x;
    I = [find(r(INQ) <= 0); EQ];
    if ninq > 0,
        J(INQ,iter) = r(INQ) <= 0;
    end  %  if ninq > 0,
    if neq > 0,
       J(EQ,iter) = ones(neq,1);
    end  %  if neq > 0,
    z = max(-r(INQ),zeros(ninq,1));
    z = [z; -r(EQ)];
% Function value:
    fval = 0.5*z'*z;
% Gradient value:
    temp = A'*z;
    gradient = temp'*temp;
    if (printlevel > 0)
       disp(' ');
       disp('---------------------------------------------------------------')
       disp(['Before iteration ',num2str(iter)]);
       disp(['            Initial function value is ',num2str(fval)])
       disp(['            Initial gradient value is ',num2str(gradient)])
       disp(['            Initial active indices: ',num2str(neq+sum(r(INQ)<0))])
    end % if (printlevel > 0)
%
%  Initialize iteration counter
%
    iters = 0;
%
%  Initialize permutation vector E if using QR methods:
%
    if (dirmeth > 1)
       E = [1:n]';
    end % if (dirmeth > 1)
%
%  Unpack tolerances:
%
    epsout = tolerances(1);
    epsline = tolerances(2);
%
%  Begin iterations if fval, gradient bigger than tolerance:
%
    while (iter <= maxits & fval > epsout & gradient > epsout)
%
%---------------------------------------------------------
% Search direction and product with A:
%---------------------------------------------------------
%
%  dirmeth == 0  means complete orthogonal decomposition each time
%  dirmeth == 1  means do singular value decomposition each time
%  dirmeth == 2  means do QR with pivoting each time
%  dirmeth == 3  means do up/downdate if possible
%
	tic;
    if (dirmeth == 0)
         [d,R,error,E] = searchdir(A,I,J(:,iter),r,0,E);
         if (printlevel > 1)
            disp(['                SD info: search direction found by COF'])
         end  % if (printlevel > 1)
    elseif (dirmeth == 1)
         [d,R,error,E] = searchdir(A,I,J(:,iter),r,1,E);
         if (printlevel > 1)
            disp(['                SD info: search direction found by SVD'])
         end  % if (printlevel > 1)
    elseif (dirmeth == 2)
         [d,R,error,E] = searchdir(A,I,J(:,iter),r,2,E);
         if (printlevel > 1)
            disp(['                SD info: search direction found by QRE'])
         end  % if (printlevel > 1)
    elseif (dirmeth == 3)
       if (iter == 1 )
           [d,R,error,E] = searchdir(A,I,J(:,iter),r,2,E);
           if (printlevel > 1)
              disp(['                SD info: search direction found by QRE'])
           end  % if (printlevel > 1)
       else % dirmeth = 3 and iter > 1 in this case
           if (m <= n ),
               [d,R,error,E] = searchdir(A,I,J(:,iter-1),r,2,E);
               if (printlevel > 1)
                  disp(['                SD info: search direction found by QRE'])
               end  % if (printlevel > 1)
           else, % try update/downdate
               if (printlevel > 1)
                  disp(['                SD info: attempting up/downdate '])
               end  % if (printlevel > 1)
               [d,R,error,E] = searchdir(A,I,J(:,iter-1),r,3,E);
               if (error == 1)
                  [d,R,error,E] = searchdir(A,I,J(:,iter-1),r,2,E);
                  if (printlevel > 1)
                     disp(['                SD info: search direction found by QRE'])
                  end  % if (printlevel > 1)
               else
                  if (printlevel > 1)
                     disp(['                SD info: search direction found by up/downdate'])
                  end  % if (printlevel > 1)
               end % if (error == 1)
           end % if (m <= n )
       end % if (iter == 1 )
    end % if (dirmeth == 1)
	timereqd = toc;

    if (printlevel > 1)
        disp(['                SD info: time required was ',num2str(timereqd)])
    end  % if (printlevel > 1)
%
% D stores all the normalized search directions.  This
% can be used to determine later if some form of conjugacy
% exists among the directions.
%
%    D = [D,d/norm(d)];
%
%
%---------------------------------------------------------
% Line search:
%---------------------------------------------------------
%
     tic;
     [stepsize,fval,linederiv,err] = linesrch(A,x,d,r,fval,epsline,linemeth,neq,printlevel);
     if (err ~= 0),
        if (printlevel > 0)
           disp(['                LS info: Error in line search; err = ',num2str(err)])
        end  % if (printlevel > 0)
        return
     end % if (err ~= 0),
	 timereqd = toc;
     if (printlevel > 1)
        disp(['                LS info: time required was ',num2str(timereqd)])
     end  % if (printlevel > 1)
% 
% Update x to x + stepsize*d;
%
     x = x + stepsize*d;
% 
% Compute new quantities at x:
%
     iter = iter + 1;
     r = b - A*x;
     I = [find(r(INQ) <= 0); EQ];
     if ninq > 0,
         J(INQ,iter) = r(INQ) <= 0;
     end  %  if ninq > 0,
     if neq > 0,
        J(EQ,iter) = ones(neq,1);
     end  %  if neq > 0,
     z = max(-r(INQ),zeros(ninq,1));
     z = [z; -r(EQ)];
     temp = A'*z;
     gradient = temp'*temp;
     if (printlevel > 0)
         disp('---------------------------------------------------------------')
         disp(['Before iteration ',num2str(iter)]);
         disp(['            function value is ',num2str(fval)])
         disp(['            gradient value is ',num2str(gradient)])
         disp(['            active indices: ',num2str(neq+sum(r(INQ)<0))])
     end % if (printlevel > 0)
%
  end % while (iter <= maxits & fval > epsout & gradient > epsout)
  iters = iter-1;
