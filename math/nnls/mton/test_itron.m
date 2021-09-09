%
% TEST_ITRON is a script that demonstrates the use of the box-constraint
%            optimization routine ITRON
%
% see also: itron.m, mtron.m, plap.m, readme.txt
%

% Number of nodes for the finite difference discretization of the
% p-Laplacian functional
N = 101;
% Sobolev index for the p-Laplacian
p = 8;

% generate the initial condition
x = linspace(0, 1, N)';
u0 = 0.3 * sin(pi*x);


% the bound constraints. On interior nodes, they are (essentially)
% -\infty \leq u_j \leq \infty. on the boundary nodes the constraint
% can be used to define the boundary condition u_0 = u_N = 0
ul = -1e300 * ones(N, 1); uu = -ul;
ul(1) = 0.0; uu(1) = 0.0;
ul(N) = 0.0; uu(N) = 0.0;

% plot initial condition
plot(x, u0, 'b-', 'Linewidth', 2.0); drawnow; hold on;

% call itron with parameters
%  fmin = -1e30 (-\infty), gtol = 1e-9
[u, F, exitflag] = itron(@(uu)(plap(uu, p)), ...
                         u0, ul, uu, -1e30, 0.1, 1e-8, 1e3);

% plot solution
plot(x, u, 'r-', 'Linewidth', 2.0); hold off;

