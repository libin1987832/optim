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
x0=BQP.x;
A=BQP.H;
b=BQP.c;
BQP.H = A'*A;
BQP.c = A'*b;
ub=BQP.ub;
tic
[BQP,output] = ipdas(BQP,option);
toc
options = optimoptions('lsqlin','Algorithm','interior-point','Display','iter');
options.Display = 'off';
options.OptimalityTolerance = 1e-15;
options.ConstraintTolerance = 1e-20;
options.MaxIterations = 200;
[m,n] =size(BQP.H);

tic;
[x1,f1,residual,exitflag,output,ff] = lsqlin(A,b,...
    [],[],[],[],[],ub,x0,options);
toc;
0.5*x1'*A'*A*x1 - (A*b)'*x1
end
