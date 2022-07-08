clc
clear
str = 'bcsstk02.mtx bcsstk04.mtx';
strindex = [1 12; 14 25];
iterstr = [98 422;984 973];
for i = 1:size(strindex,1) 
[B,rows,cols,entries,rep,field,symm] = mmread(str(1,strindex(i,1):strindex(i,2)));
n=cols;
d = sprand(n,1,1);
x1=sparse(ones(n,1)./n);
t=clock;
a = (d' * x1) / (x1' * B *x1);
[x, iter, Abpp,error] =  FLQP(a, B, d, 0, n, 2, 100, 20, 1e-6, 1e-9,1);
tf=etime(clock,t);
xBx = x' * B * x;
dx = d'* x;
ay = dx / xBx;
fa = dx - ay * xBx;
disp(['flqp:it=' num2str(iter) ',a=' num2str(ay) ',fa='...
    num2str(fa) ',abpp=' num2str(Abpp) ',nLsyst=' num2str(Abpp * iter) ',Cpu=' num2str(tf)])
figure;
semilogy(1:iter,error(1:iter));
% funxAbx = @(x) -(x'* d)/(x' * B * x);
% [xf,yval]=fmincon(funxAbx,x1,[],[],ones(1,n),1,zeros(n,1),[],[]);
% dxf = d'* xf;
% xfBxf = xf' * B * xf;
% ayf = dxf / xfBxf
end


% options.lb = zeros(1,n);  % Lower bound on the variables.
%   options.ub = [];  % Upper bound on the variables.
% 
%   % The callback functions.
%   funcs.objective        = @objective;
%   funcs.gradient         = @gradient;
%   funcs.hessian          = @hessian;
%   funcs.hessianstructure = @hessianstructure;
%   funcs.iterfunc         = @callback;
% 
%   % Set the IPOPT options.
%   options.ipopt.mu_strategy = 'adaptive';
%   options.ipopt.print_level = 0;
%   options.ipopt.tol         = 1e-7;
%   options.ipopt.max_iter    = 100;
% 
%   % Run IPOPT.
%   [x info] = ipopt(x0,funcs,options);
% 
% % ----------------------------------------------------------------------
% function f = objective (x)
%   f = 100*(x(2)-x(1)^2)^2 + (1-x(1))^2 + 90*(x(4)-x(3)^2)^2 + (1-x(3))^2 + ...
%       10.1*(x(2)-1)^2 + 10.1*(x(4)-1)^2 + 19.8*(x(2)-1)*(x(4)-1);
% 
% % ----------------------------------------------------------------------
% function g = gradient (x)
%   g(1) = -400*x(1)*(x(2)-x(1)^2) - 2*(1-x(1));
%   g(2) = 200*(x(2)-x(1)^2) + 20.2*(x(2)-1) + 19.8*(x(4)-1);
%   g(3) = -360*x(3)*(x(4)-x(3)^2) -2*(1-x(3));
%   g(4) = 180*(x(4)-x(3)^2) + 20.2*(x(4)-1) + 19.8*(x(2)-1);
%   
% % ----------------------------------------------------------------------
% function H = hessianstructure()
%   H = sparse([ 1  0  0  0 
%                1  1  0  0
%                0  0  1  0
%                0  1  1  1 ]);
% 
% % ----------------------------------------------------------------------
% function H = hessian (x, sigma, lambda)
%   H = [ 1200*x(1)^2-400*x(2)+2  0       0                          0
%         -400*x(1)               220.2   0                          0
%          0                      0       1080*x(3)^2- 360*x(4) + 2  0
%          0                      19.8   -360*x(3)                   200.2 ];
%   H = sparse(sigma*H);
%   
% % ----------------------------------------------------------------------
% function b = callback (t, f, x)
%   fprintf('%3d  %0.3g \n',t,f);
%   b = true;