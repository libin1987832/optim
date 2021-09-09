%
% PLAP implements a finite difference discretizaion of the p-Laplacian
%      functional
%
% Let I(u) = \int_0^1 |u|^p - u  dx
% Then PLAP implements the following finite difference functional:
%   I_h(u) = \sum_{j = 1}^N N^{-1} |N(u_j-u_{j-1})|^p - \sum_{j = 0}^N u_j.
% 
% Calling convention:
%    [F, G, H] = plap(u, p)
%
% Input Parameters:
%    u : N+1 column vector
%    p : Sobolev index, 1 < p < \infty (though p < 2 may be problematic)
% Output:
%    F : I_h(u)
%    G : gradient
%    H : sparse hessian
%

function [F, G, H] = plap(u, p)

N = length(u)-1;     % number of nodes
Np = N^p;            % 
d = diff(u);         % u_j - u_{j-1} (~ the gradient)

% F = \sum_{j = 1}^N  N^{-1}  | (u_j - u_{j-1}) / N |^p 
%     - \sum_{j = 0}^N N^{-1} u_j
F = (Np * sum(abs(d).^p) - sum(u)) / N; 

% compute gradient only if requested by user
if nargout > 1
    G = -ones(N+1, 1);
    G(2:(N+1)) = G(2:(N+1)) + p * Np * abs(d).^(p-2) .* d;
    G(1:N) = G(1:N) - p * Np * abs(d).^(p-2) .* d;
    G = G / N;
end

% compute hessian only if requested by user
if nargout > 2
    a = p * (p-1) * Np * abs(d).^(p-2) / N;
    H = spdiags([ [-a;0], [[a; 0]+[0; a]], [0; -a] ], [-1,0,1], N+1, N+1);
end
