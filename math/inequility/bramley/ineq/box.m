    % ---------> Usual header stuff for code development <----------
   format compact
   format short e
   dbstop if error
	warning off;
% ---------> End usual header stuff for code development <----------
%


%
%  Search direction method to use:
%
   dirmeth = 2;
   linemeth = 1;
%
%  Tolerances and problem parameters:
%
      epsout = 1.0E-12;
      epsline = 1.0E-15;
      tolerances = [epsout;epsline];
      printlevel = 2;

      n = 2;
      m = 4;
%
%  Loop over problems:
%
   for k1 = 1:1
%
%  Set up the sample problem 
%
      A = [-1 0; 1 0; 0 -1; 0 1];
      b = -ones(4,1);
      x = [10; 10];
      maxits = 2*max(m,n);
      neq = 0;
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
      disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>');
      disp(' ');
   end; % for k1 =
