function prob = random_bqp(data)
% -------------------------------------------------------------------------
% Author : Frank E. Curtis (adapted by Zheng Han)
% Date   : April 23, 2012
% Purpose: Generates the random quadratic program
%
%          (QP)  minimize   obj = 0.5*x'*H*x + x'*c + const
%               subject to  x <= xu
%
%          given the number of variables n
%          Hessian density dH.
% -------------------------------------------------------------------------
% INPUT:
% ------
%   data structure holding problem parameters desired.
%        n   number of variables.
%        dH  density of H.
%
% OUTPUT:
% -------
%   prob structure holding random problem generated.
%        H       objective quadratic term.
%        c       objective linear term.
%        const   objective constant term.
%        xu      variable upper bounds.
%        x       initial primal solution.
%        zu      initial multipliers for x <= xu.
%        obj     objective value at x.
%        Hx      product of H with initial x.
%        g       objective gradient at x, i.e. g = Hx + c.
% -------------------------------------------------------------------------

% Load necessary functions
addpath ../src/

% Set local variables
n  = data.n;
dH = data.dH;

% Generate random H(M-matrix)
H = abs(sprandsym(n,dH));
H = eye(n)*(max(eig(H))+1) - H;

% Generate the optimal x
x = randn(n,1);

% Initialize xl and xu
xu = sparse(n,1);

% Initialize zl and zu
zu = sparse(n,1);

% Set breakpoints
cl = sparse(3,1);
for i = 1:3, cl(i) = ceil(i*n/3); end;

% Approximately 1/3th of the variables will fall into these categories:
% 1: inactive, upper finite
% 2: inactive, upper infinite
% 3: active  , upper finite

% Generate upper bound data
xu(1      :cl(1)) = x(1      :cl(1)) + 5*rand(cl(1)     ,1);
xu(1+cl(1):cl(2)) = x(1+cl(1):cl(2)) + inf;
xu(1+cl(2):cl(3)) = x(1+cl(2):cl(3));

% Generate upper bound multipliers
zu(1     :cl(2)) = 0;
zu(1+cl(2):cl(3)) = rand(cl(3)-cl(2),1);

% Generate c
c = -H*x  - zu;

% Fill output structure
prob = clsBQP(H,c,xu,[],1:n);

return
