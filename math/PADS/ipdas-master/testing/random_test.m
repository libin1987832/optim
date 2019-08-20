function [BQP,output] = random_test(n,dH,varargin)
%===========================================================================
% Description	: Test wrapper for ipdas
% Parameters	: n    	      | number of variables
% 		        : dH	        | sparsity of Hessian
% 		        : method      | ssm solver options
%                           |- 'exact'
%                           |- 'bound1'
%                           |- 'bound2'(default)
%                           |- 'inexact'
%               : tol_opt 	| tolerance for kkt residual(default 1e-10)
%               : tol_iter	| maximum number of iterations(default 1e+3)
%               : tol_res_i	| tolerance for kkt residual on I set(default 1e-10)
%               : display   | display iterations or not
%                           |- 'on'(default)
%                           |- 'off'
% Author        : Zheng Han
% Date          : 01/08/2013
% Example       : option = ipdasoptimset('method','inexact','tol_iter',1000)
%===========================================================================

addpath ../src/
data.n = n;
data.dH = dH;
option  = ipdasoptimset(varargin{:});
BQP = random_bqp(data);
[BQP,output] = ipdas(BQP,option);

end
