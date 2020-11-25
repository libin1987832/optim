%
% Perform a quadratic interpolation for the minimum.
%
%
%   Randall Bramley
%   Dept of Computer Science
%   Indiana University - Bloomington
%   24 August 1993
%
	function [lambda,val,fp,info] = quadint(Ad,r,a,b,av,bv)
%
% On entry:
% --------
% a is left hand endpoint of quadratic interval
% b is right hand endpoint of quadratic interval
% av is function value at left hand endpoint
% bv is function value at right hand endpoint
%
% On return:
% ---------
% lambda is stepsize;
% val is function value;
% fp is derivative value
%
	p = 3;
	info = 0;
	s = [a; 0.5*(a+b); b];
        m = max(size(Ad));
        temp = max(s(2)*Ad - r,zeros(m,1));
        vals = [av;0.5*sum(temp'.*temp');bv];
	coeffs = [s.*s, s, ones(size(s))]\vals;
	lambda = -coeffs(2)/(2*coeffs(1));
	temp = max(lambda*Ad - r, zeros(m,1));
	val = 0.5*temp'*temp;
% As error check, if val is not smaller than that at a,
% set lambda to be midpoint of interval [a,b].
	if (val > av)
	   info = 1;
	   lambda = s(2);
	   temp = max(lambda*Ad - r, zeros(m,1));
	   val = 0.5*temp'*temp;
	end;
        fp = Ad'*temp;
%
% Another way of doing quadratic interpolation, that uses 
% known derivative values:
%     stepsize = (b*adval-a*bdval)/(adval-bdval);
