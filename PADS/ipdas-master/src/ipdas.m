function [BQP,output] = ipdas(BQP,option)
%==================================================
% Description	: Solve Bound QP of the form
%               : min 0.5x'Hx + c'x
%               : st     x <= ub
% Algorithm     : Inexact primal-dual-active-set
% Author        : Zheng Han
% Date          : 12/11/2012
%==================================================

addpath ../testing/

output.time = cputime;

% Algorithm framework
while 1

    BQP.updateVars(option.method,option.tol_res_i);
    if strcmpi(option.display,'on')
        BQP.print();
    end
    
    % Optimal solution found
    if norm(BQP.kkt,inf) < option.tol_opt
        BQP.status = 1; break;
    end
    % Excedding iteration limit
    if BQP.iter > option.tol_iter
        BQP.status = -1; break;
    end
    
    BQP.updatePartition();
end

output.time = cputime - output.time;
end