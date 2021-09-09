%
% ITRON finds a constrained minimum of a function of several variables
%
%   ITRON attempts to solve problems of the form
%        min F(x)    subject to xl <= x <= xu  (box constraints)
%
%   Calling convention (required arguments):
%
%   [x, fval, exitflag, output] = itron(fun, x0, xl, xu, fmin)
%
%   The following additional arguments are optional. Replace by [] if
%   default values should be used.
%
%   [x, fval, exitflag, output] = itron(fun, x0, xl, xu, fmin, ...
%                                 delta, gtol, maxnit, frtol, fatol, ...
%                                 cgtol, display, output_fcn)
%
%   Required Arguments:
%      fun    : function handle. Takes one argument (x) and returns
%               [F(x), G(x), H(x)], where G(x) \in R^n is the gradient
%               of F at x and H(x) is sparse and is the hessian of F at x.
%               (note that the hessian MUST be supplied)
%      x0     :  initial guess, n-column vector
%      xl, xu : lower and upper bound constraints for x. If you don't
%               have bounds, just use +- 1e300.
%      fmin   : lower bound for F(x). If TRON/mtron finds an iterate for
%               which F(x) < fmin, it terminates with a WARNING message.
%               if you don't have a lower bound, just use -1e300.
%
%   Optional Arguments:
%      delta  : initial trust region radius. 
%               [Default : 1.0]
%      gtol   : at the end of each iteration, itron, determines which
%               variables are free (= constraints are inactive) and
%               computed the inf-norm of the gradient for those variables
%               only. If the norm is less than gtol, itron terminates
%               [Default : 1e-8]
%      maxnit : Maximum number of iterations
%               [Default : 100]
%      frtol  : TRON terminates if the relative change in the objective
%               function is less than options.frtol.
%               Default : 0.0
%      fatol  : TRON terminates if the absolute change in the objective
%               function value is less than options.fatol
%               Default : 0.0
%      cgtol  : termination tolerance for PCG iteration.
%               Default : 1e-1 (as suggested by authors of TRON)
%     display : Sets the display level: 0 = no display, 1 = display only
%               after termination, 2 = at every iteration, 3 = full detail
%               Default is 2.
%  output_fcn : set to a function which can display output of the
%               optimization process at each iteration.
%
%   Output Arguments:
%           x : final iterate. This is always returned, even if ITRON quits 
%               with and error. The user should therefore always call with 
%               at least three output arguments.
%        fval : objective function value at final iterate
%    exitflag : itron was succesful if exitflag > 0, unseccesful otherwise.
%               0 : unknown error; 1 : termination as gradient tolerance
%               was reached; 2 : termination as change in objective
%               function was smaller than a tolerance; -1 : maximum number
%               of iterations reached; -2 : TRON issued a warning.
%      output : struct containing additional information about the
%               optimization process. currently only output.message is
%               implemented
%
% see also: mtron.m, test_itron.m, readme.txt
%

%
% by Christoph Ortner
% First Version: 28/11/2006
% Changes:
%   01/12/2006
%     - extended DISPLAY options
%     - print delta at each iteration (still a bug it seems)
%   19/03/2007
%     - improved computation of norm of gradient using the indfree array
%   23/04/2007
%     - trust region radius is now correctly displayed
%     - new parameters: display and output_fcn.
%   02/05/2007
%     - display number of CG iterations
%     - better formatting of output
%

%
% TODO: - xtol,
%       - first computation of normG is invalid!!!!!
%       - statistics : number of F, G, H evaluations
%

function [x, fval, exitflag, output, G, H] = ...
  itron(fun, x0, xl, xu, fmin, delta, ...
        gtol, maxnit, frtol, fatol, cgtol, ...
        display, output_fcn)


% 0 : no display
% 1 : final display
% 2 : at every iteration
% 3 : also  display every message.
DEFAULT_DISPLAY = 2;

%% plotting frequence
PLOTFREQ = 20;
plotcnt = 0;


%
% supply optional variables
%
if (nargin < 5)
    error('mtron:itron:inputmin', 'itron : At least 5 input variables required.');
end
if (nargin > 13)
    error('mtron:itron:inputmax', 'itron : At most 11 input variables.');
end    

% make all variables which were not passed empty arrays
if (nargin < 6), delta = []; end
if (nargin < 7), gtol = []; end
if (nargin < 8), maxnit = []; end
if (nargin < 9), frtol = []; end
if (nargin < 10), fatol = []; end
if (nargin < 11), cgtol = []; end
if (nargin < 12), display = []; end
if (nargin < 13), output_fcn = []; end

% give default values to all variables which are empty arrays (i.e.
% either not passed, or passed as empty arrays
if isempty(delta), delta = 1.0; end
if isempty(gtol), gtol = 1e-8; end
if isempty(maxnit), maxnit = 100; end
if isempty(frtol), frtol = 1e-12; end
if isempty(fatol), fatol =  1e-12; end
if isempty(cgtol), cgtol = 1e-1; end
if isempty(display), display = DEFAULT_DISPLAY; end

% set the actual display variable
DISPLAY = display;
% set iteration counter to 0.
nit = 0;
% set exitflag to 0
exitflag = 0;
% set the counter for PCG iterations.
ncgit = 0;

% evaluate functional at initial condition and output interation info
[F, G, H] = fun(x0);
normG = norm(G, inf); % compute_norm(G);
if (DISPLAY >= 2)
  output = sprintf('    n          F            ||G||           delta       #PCG');  
  disp(output);
  output = sprintf('--------------------------------------------------------------');  
  disp(output);
  output = sprintf('%6d    %e   %e   %e   %4d', nit, F, normG, delta, 0);
  disp(output);
end

% initialize MTRON
[task, x, info] = mtron('INIT', x0, xl, xu, fmin, ...
                        F, G, tril(H, -1), full(diag(H)), ...
                        delta, frtol, fatol, cgtol);

% start actual loop. loop as long as the task returned by mtron
% is neither 'CONVERGENCE' nor 'WARNING'.
while ( (task(1) ~= 'C') && (task(1) ~= 'W') )
    
    if (DISPLAY >= 3)
        disp(['itron-task: ', task]);
    end
    
    % if the task is 'F', it means we have to evaluate F(u) and return it
    if (task(1) == 'F');
        F = fun(x);
        [task, x, info] = mtron('STEP', F);
        delta = info.delta;
        
	% if the task is 'GH' we have to evaluate the gradient and the
    % hessian and return them.
    elseif (task(1) == 'G')
        [F, G, H] = fun(x);
        [task, x, info] = mtron('STEP', G, tril(H, -1), full(diag(H)));
        delta = info.delta;
        
    % if the task is 'NEWX', then mtron just returned a new iterate,
    % i.e. x changed
    elseif (task(1) == 'N')
        % increase iteration counter
        nit = nit + 1;
        % evaluate F, G
        [F, G] = fun(x);
        % compute norm of gradient
        normG = compute_norm(G);
        % output iteration info
        if (DISPLAY >= 2)
          output = sprintf('%6d    %e   %e   %e   %4d', nit, F, normG, ...
                           info.delta, info.ncgit - ncgit);
          disp(output);
        end
        % update number of pCG iterations
        ncgit = info.ncgit;
        
        % if the norm of G is less than a given tolerance
        if (normG < gtol)
            exitflag = 1;
            break;
        end
        
        % pass the data to an external routine to plot the progress
        if ~isempty(output_fcn)
          plotcnt = plotcnt + 1;
          if (plotcnt >= PLOTFREQ)
            output_fcn(x, F);
            plotcnt = 0;
          end
        end

        % call TRON again. Since the current task is 'NEWX', no new
        % information has to be passed on.
        [task, x, info] = mtron('STEP');
        delta = info.delta;
    end
    
    % check iteration count
    if (nit >= maxnit)
        exitflag = -1;
        break;
    end

end

% cleanup the memory used by mtron.
mtron('CLEAN');

% delete the output variable which we now need to use for
% something else.
clear output

% write F, G and H output.
fval = F;
if (nargout >= 6)
    [F, G, H] = fun(x);
end

% output termination reason.
% 1. ||G|| small
if (exitflag == 1)
    if (DISPLAY >= 1)
        disp('itron: Termination due to ||G|| less than tolerance.');
    end
    output.message = '||G||_inf less than tolerance.';
% 2. maxnit reached
elseif (exitflag == -1)
    if (DISPLAY >= 1)
        disp('itron: Maximum number of iterations reached.');
    end
    output.message = 'Maximum number of iterations reached.';
% 3. broke loop through task message:
else
    % 3.1  change in function value less than tol.
    if (task(1) == 'C')
        if (DISPLAY >= 1)
            disp('itron: Change in function value less then tol.');
            disp(task);
        end
        exitflag = 2;
        output.message = task;
	% 3.2 warning message. 
    elseif (task(1) == 'W')
        if (DISPLAY >= 1)
            disp('itron: Termination of TRON due to error.');
            disp(task);
        end
        exitflag = -2;
        output.message = task;
    end
end



%
% computes the norm over all indices where
% x does not satisfy one of the bound constraints.
%
function normG = compute_norm(G)

% get a list of indices which are free.
indfree = mtron('FREE');
indfree = indfree(find(indfree));
% compute norm over those indices only.
normG = norm(G(indfree), inf);



