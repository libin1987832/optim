%
% MTRON is a routine that interfaces the DTRON routine of the
%       TRON bound-constriant optimization package. See readme.txt for 
%       further details.
%
% There are three possible calling conventions:
%
% function [task, x, info] = mtron('INIT', x0, xl, xu, fmin, ...
%                                  F, G, H, Hdiag, ...
%                                  delta, frtol, fatol, cgtol)
% function [task, x, info] = mtron('STEP'{, ...})
% function mtron('CLEAN')
%
% 'INIT' mode:
%    Note: If mtron('INIT') is called, but mtron is already initialized
%          then 'CLEAN' is automatically performed!
%    Input Parameters:
%       x0 : initial guess
%       xl, xu : lower and upper bounds for x
%       fmin : lower bound for objective function F
%       F, G, H, Hdiag : objective function value, gradient, hessian at x0.
%       delta : initial trust region radius
%       frtol, fatol : relative and absolute termination tolerance for 
%                      change of objective function F
%       cgtol : tolerance for termination of PCG iterations
%    Output Paramters:
%       task : 'F', 'GH', 'CONV', 'NEWX', or 'WARN'
%       x : current iterate; is always passed back whether it has
%           changed or not
%       info : contains fields  'delta' and 'ncgit', where
%              delta is the current trust region radius and ncgit is the 
%              total number of PCG iterations done so far.
%
% 'STEP' mode:
%    The input parameters depend on the current value of the 'task'
%    variable which was returned by the previous call to mtron.
%       task = 'F'    : [task, x, info] = mtron('STEP', F)
%       task = 'GH'   : [task, x, info] = mtron('STEP', G, H, Hdiag)
%       task = 'NEWX' : [task, x, info] = mtron('STEP')
%       task = 'CONV' : any further call to mtron('STEP') will result in 
%                       an error.
%       task = 'WARN' : TRON issued a warning. Iteration must be stopped.
%                       any further call to mtron('STEP') will result in
%                       an error.
%    See 'INIT' mode for the values of the output parameters.
%
% 'CLEAN' mode:
%    If mtron('CLEAN') is called then all working arrays are cleaned up
%    and the function becomes ready for a new call of mtron('INIT', ...)
%    If mtron('INIT') is called when mtron is already initialized, then
%    mtron('CLEAN') is performed automatically. This is only a convenience
%    though and the user should always call mtron('CLEAN') at the end of
%    the iteration to clean up the memory.
%
% 'FREE' mode:
%    indfree = mtron('FREE'); returns an n x 1 double array indfree
%    which contains the indices of those variables which are free and
%    is otherwise filled up with zeros.
%
%
% see also: itron.m, test_mtron.m, make_mtron.m, readme.txt
%

%
% by Christoph Ortner
% First Version: 28 November 2006
% Changes:
%  01/12/2006
%    - changed 'INIT' behaviour to do automatic 'CLEAN' first
%    - returns delta in the info
%  19/03/2007
%    - added the 'FREE' option.
%  01/05/2007
%    - info is now a struct containing delta and ncgit.
%

%
% TODO: 
% - allow user to set max.num. cg iterations
% - allow user to change the trust region radius
%

