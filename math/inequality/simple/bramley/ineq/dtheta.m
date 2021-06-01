%
%     function [val] = dtheta(lambda,Ad,r,neq)
%
% Evaluate the derivative of the line search function theta at
% the point lambda.  It is assumed that Ad = A*d and r = b - A*x
% are is initialized.  Furthermore, the last me rows of A correspond
% to equality conditions in the least squares problem.
%
%   Randall Bramley
%   Dept of Computer Science
%   Indiana University - Bloomington
%   24 August 1993
%
%
      function [val] = dtheta(lambda,Ad,r,neq)
      m = size(Ad,1);
      ninq = m - neq;
%
      z1 = max(lambda*Ad(1:ninq) - r(1:ninq),zeros(ninq,1));
      z2 = lambda*Ad(ninq+1:m) - r(ninq+1:m);
      val = Ad'*[z1;z2];
      end
