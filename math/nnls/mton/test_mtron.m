
%
% Test Routine for the mtron gateway function.
%
% F(u) is the finite difference discretization of the p-Laplacian
%           \int_0^1  |u|^p - u  dx
% test_mtron minimizes F(u) subject to the constraints u(0)=u(1)=0
%
% see also: mtron.m
%

function test_mtron

% Number of nodes for the finite difference discretization of the
% p-Laplacian functional
N = 31;
% sobolev index for p-Laplacian
p = 6;

% TRON parameters
delta = 0.1;            % initial trust region radius
frtol = 0.0;            % relative and absolute tolerance for termination 
fatol = 0.0;            % due to change of F
cgtol = 1e-1;           % tolerance for the CG iteration. This value
                        % is suggested by the authors of TRON
fmin = -1e30;           % lower bound for F (usually it is best to make
                        % it a trivial lower bound
                        
gtol = 1e-8;            % gtol is a tolerance setting not used by TRON
                        % but in this routine. The iteration terminates
                        % if the norm of the gradient is less then 
                        % gtol.

% generate the initial condition
x = linspace(0, 1, N)';
u0 = sin(pi*x);

% the bound constraints. On interior nodes, they are (essentially)
% -\infty \leq u_j \leq \infty. on the boundary nodes the constraint
% can be used to define the boundary condition u_0 = u_N = 0
ul = -1e30 * ones(N, 1);
uu = -ul;
ul(1) = 0.0; uu(1) = 0.0;
ul(N) = 0.0; uu(N) = 0.0;

% plot initial condition
plot(x, u0); hold on;

% define the objective function
% note that plap implements the entire (sparse) hessian matrix, even
% though mtron requires the striclty lower triangular part and the
% diagonal separately. This is done to demonstrate how to convert
% efficiently.
ofun = @(uu)(plap(uu, p));

% evaluate functional at initial condition and output interation info
[F, G, H] = ofun(u0);
normG = norm(G(2:(N-1)), inf);
disp([' F = ', num2str(F), '    ||G|| = ', num2str(normG)]);

% initialize MTRON
% the hessian is always passed in two parts, the strictly lower
% triangular part is passed as a sparse matrix and the diagonal
% as a full vector.
[task, u, info] = mtron('INIT', u0, ul, uu, fmin, ...
                        F, G, tril(H, -1), full(diag(H)), ...
                        delta, frtol, fatol, cgtol);

% initialize iteration counter
nit = 0;

% start actual loop. loop as long as the task returned by mtron
% is neither 'CONVERGENCE' nor 'WARNING'.
while ( (task(1) ~= 'C') && (task(1) ~= 'W') )
    
    % if the task is 'F', it means we have to evaluate F(u) and return it
    if (task(1) == 'F');
        F = ofun(u);
        [task, u, info] = mtron('STEP', F);
        
	% if the task is 'GH' we have to evaluate the gradient and the
    % hessian and return them.
    elseif (task(1) == 'G')
        [F0, G, H] = ofun(u);
        [task, u, info] = mtron('STEP', G, tril(H, -1), full(diag(H)));
        
    % if the task is 'NEWX', we plot the current iterate, evaluate the
    % gradient and the function and display the iteration info.
    elseif (task(1) == 'N')
        nit = nit + 1;
        plot(x, u);
        [F, G] = ofun(u);
        % (don't take first and last entry in norm since those are
        %  restricted by the bound constraints.)
        normG = norm(G(2:(N-1)), inf);
        disp([num2str(nit), ': F = ', num2str(F), ...
            '    ||G|| = ', num2str(normG)]);
        title([' ||G|| = ', num2str(normG)]);
        drawnow;
        
        % if the norm of G is less than a given tolerance
        if (normG < 1e-6)
            break;
        end

        % call TRON again. Since the current task is 'NEWX', no new
        % information has to be passed on.
        [task, u, info] = mtron('STEP');
    end
    
end

% cleanup the memory used by mtron.
mtron('CLEAN');

% output termination reason.
if (task(1) == 'C')
    disp('Termination due to change in function value less then tol.');
    disp(task);
elseif (task(1) == 'W')
    disp('Termination of TRON due to error.');
    disp(task);
else
    disp('Termination due to ||G|| less then tol');
end

% plot final solution
plot(x, u, 'r-', 'Linewidth', 2.0);
hold off;
