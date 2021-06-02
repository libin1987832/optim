%
%  function [stepsize,newfval,newderiv,err] = linesrch(A,x,d,r,fval,epsline,linemeth,neq);
%
%   Line search for Han's method for linear inequalities.
%
%   Finds lambda minimizing theta(lambda) = || (A(x + lambda*d) - b)_+ ||,
%     where the subscript + indicates gangster projection onto
%     the nonnegative orthant in m-dimensional space.
%
%   The calling sequence has:
%      A = matrix; assumed now to consist of all inequalities.
%      x = current iterate
%      d = search direction
%      r = b - A*x
%      fval = 1/2 || (Ax - b)_+ ||^2; function value at current iterate
%      epsline = tolerance to use in binary search methods.
%      linemeth = method to use; this has values
%            => 1  Search knot points until quadratic interval is found, then interpolate
%            => 2  Perform method 1, but if derivative of theta is bigger than epsline,
%                        refine the stepsize using binary search.
%            => 3  Use brute force binary search throughout.
%      neq = number of equality rows in A; they are assumed to be the last neq rows, with
%            0 .le. neq .le. m
%
%   On return:
%      stepsize is what you think it is
%      newfval is function value at stepsize
%      newderiv is line search function derivative value at stepsize
%      err  = 0 if all's OK
%           = 1 if d is not a search direction
%
%   Randall Bramley
%   Dept of Computer Science
%   Indiana University - Bloomington
%   24 August 1993
%   Revised 1 October 1993 and countless other times.
%
%   Revised Mon Dec 23 10:01:01 EST 2002
%		to remove flop counts for Matlab 6.X
%
%---------------------------------------------------------
%
  function [stepsize,newfval,newderiv,err] = ...
            linesrch(A,x,d,r,fval,epsline,linemeth,neq,printlevel);

  [m,n] = size(A);
% Number of inequality rows of A:
  ninq = m - neq;
% Index sets corr to inequality and equality rows of A:
  INQ = 1:ninq;
  EQ = ninq+1:m;
  INQ = INQ';
  EQ = EQ';
%
  Ad = A*d;
  err = 0;
  z = max(r(INQ),zeros(ninq,1));
  z = [z;r(EQ)];
%
% 
% Initialize search interval left hand endpoint to be a = 0:
% 
  anode = 0;
% Function value at lefthand endpoint a:
  aval = fval;
% Derivative value at lefthand endpoint a:
  adval = -Ad'*z;
% Check to make sure d is a descent direction:
  if (adval > 0) 
      disp('                LS info: Error in line search; theta increasing at 0')
      err = 1;
      return;
  end;
%
% Handle cases needing quadratic interpolation:
%
  if (linemeth == 1 | linemeth == 2)
% 
%  Find knot points which are positive and well-defined:
% 
      knots = r(INQ)./Ad(INQ);
      ind = find(~(isnan(knots)) & ~(isinf(knots)) &  knots>0);
      knots = knots(ind);
      numknots = max(size(knots));
      if (numknots <= 0)
% No positive knot points; use [0,2] as search interval
          bnode = 2;
% Find value at right hand endpoint b:
          z = max(-bnode*Ad(INQ)+r(INQ),zeros(ninq,1));
          z = [z;bnode*Ad(EQ)-r(EQ)];
%
          bval = 0.5*z'*z;
% Perform quadratic interpolation to get minimum
          [stepsize,newfval,newderiv,info] = quadint(Ad,r,anode,bnode,aval,bval);
          if (info ~= 0 & printlevel > 0),
              disp(['     LS info: *** In quadint, interpolated value set to midpoint of'])
              disp(['     LS info: ***           search interval since t(a) > t(lambda) '])
          end % if (info ~= 0),
      else
% Else, there were positive knot points ....
%
% Positive knot points exist; search through them to find unimodal interval:
% Start with largest knot point, given by ind:
%
          knots = sort(knots);
          [bnode,ind] = max(knots);
          z = max(knots(ind)*Ad(INQ)-r(INQ),zeros(ninq,1));
          z = [z;knots(ind)*Ad(EQ)-r(EQ)];
          bval = 0.5*z'*z;     %  bdval = Ad'*z;
          for i = 1:numknots;
% Only need to check knot if it is smaller than current b.
             if (knots(i) < bnode )
                 z = max(knots(i)*Ad(INQ)-r(INQ),zeros(ninq,1));
                 z = [z;knots(i)*Ad(EQ)-r(EQ)];
                 temp = Ad'*z;  % temp is the derivative at current knot point.
                 if (temp >= 0 & knots(i) < bnode)
                    if (i>1)
                        anode = knots(i-1);
                        za = max(anode*Ad(INQ)-r(INQ),zeros(ninq,1));
                        za = [za;anode*Ad(EQ)-r(EQ)];
                        aval = 0.5 * za' * za;
                    end
                    bnode = knots(i);
                    bval = 0.5*z'*z;      %   bdval = temp;
                 end; % if
             end;  % if
          end;  % for i
      end % if (numknots <= 0)
%
% Guard against humongous bnode values by arbitrarily restricting bnode:
%
          if (bnode > 100)
              if (printlevel > 0)
              	disp(['                LS info: *** bnode in linesrch is ',num2str(bnode)])
              	disp(['                LS info: ***       bnode is being restricted to 10'])
			  end;
              bnode = 10;
              z = max(bnode*Ad(INQ)-r(INQ),zeros(ninq,1));
              z = [z;bnode*Ad(EQ)-r(EQ)];
              bval = 0.5*z'*z;
          end;  % if
% Now [anode,bnode] is an interval on which theta is quadratic; perform interpolation
          [stepsize,newfval,newderiv,info] = quadint(Ad,r,anode,bnode,aval,bval);
          if (info ~= 0 & printlevel > 0),
              disp(['                LS info: *** In quadint, interpolated value set to midpoint of'])
              disp(['                LS info: ***           search interval since t(a) > t(lambda) '])
          end % if (info ~= 0),
  end  % case where methods 1 and 2 perform quadratic interpolation approach.
%
% Now if linemeth = 3, perform binary search to get interval, and then
% binary search to find minimum in that interval.
% If linemeth  = 2, check derivative and do second step if it is too large:
%
  if (linemeth == 3) 
%
% Initialize search interval to be (0,2].
%
      bnode = 2;
      dprev = max(1,-adval);
      bdval = dtheta(bnode,Ad,r,neq);
%
% If derivative of line search function is negative at bnode,
% must increase value of bnode to where derivative is nonnegative
% to be sure to include minimum point in search interval:
%
      if bdval < 0;
         kk = 1;
         while bdval < 0 & kk < 8;  % double size of b until reach 2^8 if necessary:
          if (printlevel > 0),
             disp(['                LS info: *** Line search function is descending at lambda = ',num2str(bnode)])
             disp(['                LS info: ***    Derivative of theta(bnode) is ',num2str(bdval)])
             disp('                LS info: ***    Upper bound on interval is now doubled')
%             disp('                 LS info: ***    Hit return to continue')
%             pause;
			end % if printlevel > 0
             kk = kk + 1;
             bnode = 2*bnode;
             bdval = dtheta(bnode,Ad,r,neq);
         end; % while bdval < 0 & kk < 8;
      end;  % if bdval < 0
  end % if (linemeth == 3) 
%
% Completed finding interval in case of linemeth = 3.  Now for
% linemeth = 2 or 3, perform binary search if needed.  For linemeth = 2,
% just use the value given by stepsize if it has positive derivative.
% Otherwise, use bnode:
%
  if (linemeth == 2 | linemeth ==3) 
      if (linemeth == 2),
         if (newderiv > 0 ),
            bnode = stepsize;
            bdval = newderiv;
         else,
            bdval = dtheta(bnode,Ad,r,neq);
         end % if (temp > 0 ),
      end  % if (linemeth == 2),
%
%  First check; if derivative at bnode is small enough, just use it as
%  stepsize.  Note that this may not give smallest value of lambda 
%  possible.
%
      if bdval < epsline;
         stepsize = bnode;
         newfval = bval;
         newderiv = bdval;
      else  % Must do search inside interval:
%
% Binary search with upper limit of 200 iterations.  This
% is used as a safeguarded approach to finding the stepsize.
%
         kk = 1;
         while bnode - anode > epsline & kk < 200;
             stepsize = (anode + bnode)/2;
             newderiv = dtheta(stepsize,Ad,r,neq);
             if newderiv < 0;
                anode = stepsize;
             else
                bnode = stepsize;
             end; % if newderiv < 0
             kk = kk + 1;
         end; % while bnode - anode > epsline & kk < 300
		 if (printlevel > 0)
           disp(['                LS info: Binary search required ',num2str(kk),' trial values'])
		 end;
         newfval = theta(stepsize,Ad,r);
      end; % if bdval < epsline
%
  end  % if (linemeth == 2 | linemeth ==3) 

%--------------------------------------------------------------------------
%
% Display step information
%
	 if (printlevel > 0)
     	disp(['                LS info: Stepsize = ',num2str(stepsize)]);
     	disp(['                LS info: Search interval width = ',num2str(bnode - anode)]);
%     	disp(['                LS info: Distance moved = ',num2str(norm(d)*stepsize)]);
	 end;
%
%--------------------------------------------------------------------------
% Display in a plot the graph of the line search function:
%
%  showstep;
%  disp('Hit return to continue')
%  pause;
%
% Display in a plot the graph of the derivative of the line search function:
%
%  showderiv;
%  disp('Hit return to continue')
%  pause;
%--------------------------------------------------------------------------
  return
