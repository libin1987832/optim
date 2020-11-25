%
%     function [val] = theta(lambda,Ad,r)
%
% Evaluate the line search function theta at the point lambda.
% It is assumed that Ad = A*d and r = b - A*x is initialized.
%
%
%   Randall Bramley
%   Dept of Computer Science
%   Indiana University - Bloomington
%   24 August 1993
%
      function [val] = theta(lambda,Ad,r)
      m = size(Ad);
%
      s = max(lambda*Ad -r,zeros(m,1));
      val = 0.5*s'*s;
      end
