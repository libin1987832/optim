function [x, y] = pcnnls(A,b,x0,niter)
%[x, y] = pcnnls(A,b,x0,niter);
%solves the linear least squares problem with nonnegative variables using the predictor-corrector algorithm in [1].
%Input:
%   A:      [MxN] matrix 
%   b:      [Mx1] vector
%   x0:     [Nx1] vector of initial values. x0 > 0. Default value: ones(n,1)
%   niter:  Number of iterations. Default value: 10
%Output
%   x:      solution
%   y:      complementary solution
%
% [1] Portugal, Judice and Vicente, A comparison of block pivoting and
% interior point algorithms for linear least squares problems with
% nonnegative variables, Mathematics of Computation, 63(1994), pp. 625-643
%
%Uriel Roque
%28.07.2005
%02.05.2006
[m,n] = size(A);
if nargin < 3
    x0 = ones(n,1);
    niter = 10;
elseif nargin < 4
    niter = 10;
end
if isempty(x0)
    x0 = ones(n,1);
elseif any(x0 <= 0)
    disp('Error. The initial vector should be nonzero');
    x = [];
    y = [];
    return
end
AtA = A'*A;
Atb = A'*b;
In = eye(n);
e = ones(n,1);
%Step 0
y0 = x0;
k = 0;
TOL1 = n*eps;
TOL2 = n*sqrt(eps);
xk = x0;
yk = y0;
noready = 1;
while noready
    
    k= k+1;
    
    %Step 1
    Xk = diag(xk);
    Yk = diag(yk);
    C = [Xk Yk; AtA -In];
    d = [-Xk*Yk*e; -AtA*Xk*e+Yk*e+Atb];
    uv = C\d;
    uk(:,1) = uv(1:n);
    vk(:,1) = uv(n+1:2*n);
    
    %Step 2
    T1 = min(-xk(uk < 0) ./ uk(uk < 0)); 
    T2 = min(-yk(vk < 0) ./ vk(vk < 0));
    theta = 0.99995 * min(T1,T2);
    if isempty(theta)
        noready = 0;
        theta = 0;
    end
    
    %Step 3
    mk = (xk + theta*uk)'*(yk + theta*vk) / (n^2);
    Uk = diag(uk);
    Vk = diag(vk);
    d = [-Xk*Yk*e+mk*e-Uk*Vk*e; -AtA*Xk*e+Yk*e+Atb];
    zw = C\d;
    zk(:,1) = zw(1:n);
    wk(:,1) = zw(n+1:2*n);
    
    %Step 4
    T1 = min(-xk(zk < 0) ./ zk(zk < 0)); 
    T2 = min(-yk(wk < 0) ./ wk(wk < 0));
    theta = 0.99995 * min(T1,T2);
    if isempty(theta)
        theta = 0;
        noready = 0;
    end
    xk = xk + theta*zk;
    yk = yk + theta*wk;
    
    %Step 5
    if (xk'*yk < TOL1) & (norm(AtA*xk-Atb-yk) < TOL2)
        noready = 0;
    elseif k == niter
        noready = 0;
    end
    
end
x = xk;
y = yk;