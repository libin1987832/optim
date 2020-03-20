A=diag([10^15,1]);
% tf = issymmetric(A);

d = eig(A);
tol=length(d)*eps(max(d))
isposdef = all(d > tol)
issemidef = all(d > -tol)