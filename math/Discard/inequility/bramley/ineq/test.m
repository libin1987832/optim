% ---------> Usual header stuff for code development <----------
   format compact
   format short e
   dbstop if error
	warning off;
%   !rm stuff
%   diary stuff
% ---------> End usual header stuff for code development <----------
%

%-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
%  Test the linear inequality solver matlab code.  This driver sets
%  up random problems in a loop (on k1), but can easily be modified
%  to solve another problem by, e.g., substituting explicitly set
%  up matrices or by reading in the data from a file.  The routine
%  uses menus and keyboard input variables to select the methods to
%  use for linesearch and search direction finding procedures.  By
%  default it will create a random 40x20 matrix A and corresponding
%  right hand side vector b.
%
%   Randall Bramley and Beata Winnicka
%   Dept of Computer Science
%   Indiana University - Bloomington
%   1 October 1993 
%-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
      disp('This driver will create a test problem of order 40 x 20 and')
      disp('call the inequality least squares solver.  The user selects')
      disp('a random seed to create the problem initially and the ')
      disp('search direction and line search methods to use.')
  
%
%  Random number generator seed:
%
   seed = 6;
   seed = input('Random seed to use: ');
   randn('seed',seed);
%
%  Search direction method to use:
%
   dirmeth = menu('Search direction method:', ...
             'Complete orthogonal decomposition', ...
             'Singular value decomposition', ...
             'QR with column pivoting', ...
             'QR with update/downdate');
	% Fix a bug introduced by shift to Matlab 6.x: numbers start from 1, not 0
	dirmeth = dirmeth - 1;
%
%  Line search method to use:
%
   linemeth = menu('Line search method:', ...
             'Knot search + quadratic interpolation', ...
             'Knot search + quadratic interpolation + binary search', ...
             'Binary search');
%
%  Tolerances and problem parameters:
%
      epsout = 1.0E-12;
      epsline = 1.0E-15;
      tolerances = [epsout;epsline];
      printlevel = 2;

      n = 20;
      m = 2*n;
%
%  Loop over problems:
%
	upper_limit = 1;
	Iterations = zeros(upper_limit,1);
   for k1 = 1:upper_limit
%
%  Set up the sample problem 
%
      A = randn(m,n);
      b = randn(m,1);
      x = randn(n,1);
      maxits = 2*max(m,n);
      neq = 0;
%      neq = input('Number of equality rows in matrix: ');
%
      disp(' ');
      disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>');
      disp(['Problem number ',num2str(k1)]);
      disp(['Number of rows: ',num2str(m)]);
      disp(['Number of cols: ',num2str(n)]);
      disp(['Tolerance on gradient: ',num2str(epsout)]);
      disp(['Tolerance on line search: ',num2str(epsline)]);
      disp(['Max iterations allowed: ',num2str(maxits)]);
      disp(['Number of equality rows in A: ',num2str(neq)]);
      disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>');
      disp(' ');
%
      [x,fval,gradient,err,J,D,iters] = ineq(A,x,b,tolerances,dirmeth,linemeth,maxits,printlevel,neq);
%
	  Iterations(k1) = iters;
      disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>');
      disp(' ');
   end; % for k1 =
	hist(Iterations);
	figure;
	plot(sort(Iterations),'+');
